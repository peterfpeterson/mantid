// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidCurveFitting/Functions/GaussianComptonProfile.h"
#include <cxxtest/TestSuite.h>

#include "ComptonProfileTestHelpers.h"
#include "MantidAPI/FunctionDomain1D.h"

using Mantid::CurveFitting::Functions::ComptonProfile;
using Mantid::CurveFitting::Functions::GaussianComptonProfile;

class GaussianComptonProfileTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static GaussianComptonProfileTest *createSuite() { return new GaussianComptonProfileTest(); }
  static void destroySuite(GaussianComptonProfileTest *suite) { delete suite; }

  void test_Name_Is_As_Expected() {
    // These are used in scripts so should not change!
    Mantid::API::IFunction_sptr profile = createFunction();
    TS_ASSERT_EQUALS("GaussianComptonProfile", profile->name());
  }

  void test_Initialized_Function_Has_Expected_Parameters_In_Right_Order() {
    Mantid::API::IFunction_sptr profile = createFunction();
    static const size_t nparams(3);

    auto currentNames = profile->getParameterNames();
    const size_t nnames = currentNames.size();
    TS_ASSERT_EQUALS(nparams, nnames);

    if (nnames == nparams) {
      const char *expectedParams[nparams] = {"Mass", "Width", "Intensity"};
      for (size_t i = 0; i < nnames; ++i) {
        TS_ASSERT_EQUALS(expectedParams[i], currentNames[i]);
      }
    }
  }

  void test_Function_Has_One_Intensity_Coefficient() {
    std::shared_ptr<ComptonProfile> profile = createFunction();

    auto intensityIndices = profile->intensityParameterIndices();
    TS_ASSERT_EQUALS(1, intensityIndices.size());
  }

  void test_Expected_Results_Returned_Given_Data() {
    using namespace Mantid::API;

    auto func = createFunctionWithParamsSet();
    double x0(370.0), x1(371.0),
        dx(0.5); // chosen to give put us near the peak for this mass & spectrum
    auto testWS =
        ComptonProfileTestHelpers::createTestWorkspace(1, x0, x1, dx, ComptonProfileTestHelpers::NoiseType::None);
    auto &dataX = testWS->dataX(0);
    using std::placeholders::_1;
    std::transform(dataX.begin(), dataX.end(), dataX.begin(),
                   std::bind(std::multiplies<double>(), _1, 1e-06)); // to seconds
    func->setMatrixWorkspace(testWS, 0, dataX.front(), dataX.back());
    FunctionDomain1DView domain(dataX.data(), dataX.size());
    FunctionValues values(domain);

    TS_ASSERT_THROWS_NOTHING(func->function(domain, values));
    const double tol(1e-6);
    TS_ASSERT_DELTA(0.1048941000, values.getCalculated(0), tol);
    TS_ASSERT_DELTA(0.1044889285, values.getCalculated(1), tol);
    TS_ASSERT_DELTA(0.1029765223, values.getCalculated(2), tol);
  }

private:
  std::shared_ptr<GaussianComptonProfile> createFunctionWithParamsSet() {
    auto func = createFunction();
    func->setParameter("Mass", 30.0);
    func->setParameter("Intensity", 4.0);
    func->setParameter("Width", 13.0);
    func->setUpForFit();
    return func;
  }

  std::shared_ptr<GaussianComptonProfile> createFunction() {
    auto profile = std::make_shared<GaussianComptonProfile>();
    profile->initialize();
    return profile;
  }
};
