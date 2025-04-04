// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidGeometry/Instrument/XMLInstrumentParameter.h"
#include "MantidGeometry/IComponent.h"
#include "MantidGeometry/muParser_Silent.h"
#include "MantidKernel/Exception.h"
#include "MantidKernel/LogParser.h"
#include "MantidKernel/Logger.h"
#include "MantidKernel/TimeROI.h"
#include "MantidKernel/TimeSeriesProperty.h"
#include <boost/regex.hpp>
#include <ctime>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>

namespace Mantid::Geometry {
namespace {
Kernel::Logger g_log("XMLInstrumentParameter");
}

using namespace Kernel;

/** Constructor
 *  @param logfileID :: The logfile id -- the part of the file name which
 * identifies the log
 *  @param value :: Rather then extracting value from logfile, specify a value
 * directly
 *  @param paramName :: The name of the parameter which will be created based on
 * the log values
 *  @param type :: The type
 *  @param extractSingleValueAs :: Describes the way to extract a single value
 * from the log file( average, first number, etc)
 *  @param eq :: muParser equation to calculate the parameter value from the log
 * value
 *  @param comp :: The pointer to the instrument component
 *  @param interpolation :: The pointer to the interpolation class
 *  @param formula :: The string formula to apply
 *  @param formulaUnit :: The unit that the formula requires the input value in
 *  @param resultUnit :: The unit of the result of the formula
 *  @param tie :: What to tie the value to
 *  @param constraint :: The constraint associated with this parameter
 *  @param penaltyFactor :: The level of penalty associated with the constraint
 *  @param fitFunc :: What fit function this applies to
 *  @param angleConvertConst :: angle conversion constant?????
 *  @param description :: text description of the parameter
 *  @param visible :: whether the parameter should be visible in InstrumentViewer
 */
XMLInstrumentParameter::XMLInstrumentParameter(std::string logfileID, std::string value,
                                               std::shared_ptr<Kernel::Interpolation> interpolation,
                                               std::string formula, std::string formulaUnit, std::string resultUnit,
                                               std::string paramName, std::string type, std::string tie,
                                               std::vector<std::string> constraint, std::string &penaltyFactor,
                                               std::string fitFunc, std::string extractSingleValueAs, std::string eq,
                                               const Geometry::IComponent *comp, double angleConvertConst,
                                               const std::string &description, std::string visible)
    : m_logfileID(std::move(logfileID)), m_value(std::move(value)), m_paramName(std::move(paramName)),
      m_type(std::move(type)), m_tie(std::move(tie)), m_constraint(std::move(constraint)),
      m_penaltyFactor(penaltyFactor), m_fittingFunction(std::move(fitFunc)), m_formula(std::move(formula)),
      m_formulaUnit(std::move(formulaUnit)), m_resultUnit(std::move(resultUnit)),
      m_interpolation(std::move(interpolation)), m_extractSingleValueAs(std::move(extractSingleValueAs)),
      m_eq(std::move(eq)), m_component(comp), m_angleConvertConst(angleConvertConst), m_description(""),
      m_visible(std::move(visible)) {
  if (!description.empty()) { // remove multiple spaces
    static const boost::regex re("\\s+");
    std::string desc = boost::regex_replace(description, re, " ");
    (const_cast<std::string *>(&m_description))->assign(desc);
  }
}

/** Returns the parameter value.
 * This interprets the XML parameter specification in order to do one of these
 *things:
 *  - Calculate an equation result, if specified
 *  - Interpolate the value, if desired.
 *  - Just extract the value (perhaps the man or just the n-th position) and
 *return that.
 *
 *  @param logData :: Data in logfile
 *  @param roi :: TimeROI object to get active time
 *  @return parameter value
 *
 *  @throw InstrumentDefinitionError Thrown if issues with the content of XML
 *instrument definition file
 */
double XMLInstrumentParameter::createParamValue(const TimeSeriesProperty<double> *logData,
                                                const Kernel::TimeROI *roi) const {
  // If this parameter is a <look-up-table> or <formula> return 0.0. Such
  // parameter types are
  // associated with 'fitting' parameters. In some sense this method should
  // never be called
  // for such parameters since they return values from other attributes and only
  // during (or just
  // before) 'fitting'. But include the statement below for safety.

  if (!m_formula.empty() || m_interpolation->containData())
    return 0.0;

  // also this method should not be called when parameter is of 'string' type.
  // Display
  // an error and return 0.0

  if (m_type == "string") {
    g_log.error() << "XMLInstrumentParameter::createParamValue has been called "
                     "with a 'string' parameters.\n"
                  << "Return meaningless zere value.";
    return 0.0;
  }

  double extractedValue = 0.0;

  // Get value either as directly specified by user using the 'value' attribute
  // or through
  // a logfile as specified using the 'logfile-id' attribute. Note if both
  // specified 'logfile-id'
  // takes precedence over the 'value' attribute

  if (!m_logfileID.empty()) {
    // get value from time series

    using StatisticsMapType = std::map<std::string, Kernel::Math::StatisticType>;
    StatisticsMapType statistics_types;
    statistics_types.emplace("first_value", Kernel::Math::FirstValue);
    statistics_types.emplace("last_value", Kernel::Math::LastValue);
    statistics_types.emplace("maximum", Kernel::Math::Maximum);
    // statistics_types.emplace("mean", Kernel::Math::Mean);
    // //TODO, would conflict with the existing "mean" flag, which corresponds
    // to time_averaged_mean
    statistics_types.emplace("median", Kernel::Math::Median);
    statistics_types.emplace("minimum", Kernel::Math::Minimum);
    StatisticsMapType::const_iterator statisics_choice = statistics_types.find(m_extractSingleValueAs);
    const bool bUsingStandardStatistics = statisics_choice != statistics_types.end();

    if (m_extractSingleValueAs == "mean") {
      extractedValue = timeMean(logData);
    } else if (bUsingStandardStatistics) {
      extractedValue = logData->extractStatistic((*statisics_choice).second, roi);
    }
    // Looking for string: "position n", where n is an integer and is a 1-based
    // index
    else if (m_extractSingleValueAs.starts_with("position") && m_extractSingleValueAs.size() >= 10) {
      std::stringstream extractPosition(m_extractSingleValueAs);
      std::string dummy;
      int position;
      extractPosition >> dummy >> position;

      extractedValue = logData->nthValue(position - 1);
    } else {
      throw Kernel::Exception::InstrumentDefinitionError(
          std::string("extract-single-value-as attribute for <parameter>") + " element (eq=" + m_eq +
          ") in instrument definition file is not recognised.");
    }
  } else {
    try {
      extractedValue = boost::lexical_cast<double>(m_value);
    } catch (boost::bad_lexical_cast &) {
      throw Kernel::Exception::InstrumentDefinitionError(std::string("<parameter> with name ") + m_paramName +
                                                         " much be set to a number,\n" +
                                                         "unless it is meant to be a 'string' parameter.");
    }
  }

  // Check if m_eq is specified if yes evaluate this equation

  if (m_eq.empty())
    return extractedValue;

  size_t found;
  std::string equationStr = m_eq;
  found = equationStr.find("value");
  if (found == std::string::npos) {
    throw Kernel::Exception::InstrumentDefinitionError(
        std::string("Equation attribute for <parameter>") + " element (eq=" + m_eq +
        ") in instrument definition file must contain the string: \"value\"." +
        ". \"value\" is replaced by a value from the logfile.");
  }

  std::stringstream readDouble;
  readDouble << extractedValue;
  std::string extractedValueStr = readDouble.str();
  equationStr.replace(found, 5, extractedValueStr);

  // check if more than one 'value' in m_eq

  while (equationStr.find("value") != std::string::npos) {
    found = equationStr.find("value");
    equationStr.replace(found, 5, extractedValueStr);
  }

  try {
    mu::Parser p;
    p.SetExpr(equationStr);
    return p.Eval();
  } catch (mu::Parser::exception_type &e) {
    throw Kernel::Exception::InstrumentDefinitionError(
        std::string("Equation attribute for <parameter>") + " element (eq=" + m_eq +
        ") in instrument definition file cannot be parsed." + ". Muparser error message is: " + e.GetMsg());
  }
}

} // namespace Mantid::Geometry
