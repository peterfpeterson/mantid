// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2025 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidLegacyNexus/DllConfig.h"
#include "MantidLegacyNexus/NeXusFile_fwd.h"
#include "MantidLegacyNexus/napi.h"

#include <boost/container/vector.hpp>
#include <ios>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Mantid {
namespace LegacyNexus {

/** C++ implementation of NeXus classes.

@author Roman Tolchenov, Tessella plc
@date 28/05/2009
*/

/** Structure for keeping information about a Nexus data set,
 *  such as the dimensions and the type
 */
struct NXInfo {
  NXInfo() : nxname(), rank(0), dims(), type(NXnumtype::BAD), stat(NXstatus::NX_ERROR) {}
  std::string nxname; ///< name of the object
  int rank;           ///< number of dimensions of the data
  int dims[4];        ///< sizes along each dimension
  NXnumtype type;     ///< type of the data, e.g. NX_CHAR, NXnumtype::FLOAT32, see napi.h
  NXstatus stat;      ///< return status
  operator bool() { return stat == NXstatus::NX_OK; } ///< returns success of an operation
};

/**  Information about a Nexus class
 */
struct NXClassInfo {
  NXClassInfo() : nxname(), nxclass(), datatype(NXnumtype::BAD), stat(NXstatus::NX_ERROR) {}
  std::string nxname;  ///< name of the object
  std::string nxclass; ///< NX class of the object or "SDS" if a dataset
  NXnumtype datatype;  ///< NX data type if a dataset, e.g. NX_CHAR, NXnumtype::FLOAT32, see
  /// napi.h
  NXstatus stat;                                      ///< return status
  operator bool() { return stat == NXstatus::NX_OK; } ///< returns success of an operation
};

/**
 * LoadNexusProcessed and SaveNexusProcessed need to share some attributes, put
 * them at
 * namespace level here
 */
/// Default block size for reading and writing processed files
const int g_processed_blocksize = 8;

/**  Nexus attributes. The type of each attribute is NX_CHAR
 */
class MANTID_LEGACYNEXUS_DLL NXAttributes {
public:
  int n() const { return int(m_values.size()); }         ///< number of attributes
  std::string operator()(const std::string &name) const; ///< returns the value of attribute with name name
  void set(const std::string &name,
           const std::string &value); ///< set the attribute's value
private:
  std::map<std::string, std::string> m_values; ///< the list of attributes
};

/// Forward declaration
class NXClass;

/**  The base abstract class for NeXus classes and data sets.
 *    NX classes and data sets are defined at www.nexusformat.org
 */
class MANTID_LEGACYNEXUS_DLL NXObject {
  friend class NXDataSet; ///< a friend class declaration
  friend class NXClass;   ///< a friend class declaration
  friend class NXRoot;    ///< a friend class declaration
public:
  // Constructor
  NXObject(const NXhandle fileID, const NXClass *parent, const std::string &name);
  virtual ~NXObject() = default;
  /// Return the NX class name for a class (HDF group) or "SDS" for a data set;
  virtual std::string NX_class() const = 0;
  // True if complies with our understanding of the www.nexusformat.org
  // definition.
  // virtual bool isStandard()const = 0;
  /// Returns the absolute path to the object
  std::string const &path() const { return m_path; }
  /// Returns the name of the object
  std::string name() const;
  /// Attributes
  NXAttributes attributes;
  /// Nexus file id
  NXhandle m_fileID;

protected:
  std::string m_path; ///< Keeps the absolute path to the object
  bool m_open;        ///< Set to true if the object has been open
private:
  NXObject() : m_fileID(), m_open(false) {} ///< Private default constructor
  void getAttributes();
};

/** Abstract base class for a Nexus data set. A typical use include:
 *  <ul>
 *       <li>Creating a dataset object using either the concrete type
 * constructor or specialized methods of NXClass'es</li> <li>Opening the dataset
 * with open() method. Specialized NXClass creation methods call open()
 * internally (so no need to call it again).</li> <li>Loading the data using
 * load(...) method. The data can be loaded either in full or by chunks of
 * smaller rank (dimension)</li>
 *  </ul>
 *  There is no need to free the memory allocated by the NXDataSet as it is done
 * at the destruction.
 */
class MANTID_LEGACYNEXUS_DLL NXDataSet : public NXObject {
public:
  // Constructor
  NXDataSet(const NXClass &parent, const std::string &name);
  /// NX class name. Returns "SDS"
  std::string NX_class() const override { return "SDS"; }
  /// Opens the data set. Does not read in any data. Call load(...) to load the
  /// data
  void open();
  /// Opens datasets faster but the parent group must be already open
  void openLocal();
  /// Returns the rank (number of dimensions) of the data. The maximum is 4
  int rank() const { return m_info.rank; }
  /// Returns the number of elements along i-th dimension
  int dims(int i) const { return i < 4 ? m_info.dims[i] : 0; }
  /// Returns the number of elements along the first dimension
  int dim0() const;
  /// Returns the number of elements along the second dimension
  int dim1() const;
  /// Returns the number of elements along the third dimension
  int dim2() const;
  /// Returns the number of elements along the fourth dimension
  int dim3() const;
  /// Returns the name of the data set
  std::string name() const { return m_info.nxname; } // cppcheck-suppress returnByReference
  /// Returns the Nexus type of the data. The types are defied in napi.h
  NXnumtype type() const { return m_info.type; }
  /** Load the data from the file. Calling this method with all default
   * arguments makes it read in all the data.
   */
  virtual void load() {};

protected:
  void getData(void *data);

private:
  NXInfo m_info; ///< Holds the data info
};

template <typename T>
using container_T = std::conditional_t<std::is_same<T, bool>{}, boost::container::vector<bool>, std::vector<T>>;

/**  Templated class implementation of NXDataSet. After loading the data it can
 * be accessed via operators () and [].
 */
template <class T> class NXDataSetTyped : public NXDataSet {

public:
  /**  Constructor.
   *   @param parent :: The parent Nexus class. In terms of HDF it is the group
   * containing the dataset.
   *   @param name :: The name of the dataset relative to its parent
   */
  NXDataSetTyped(const NXClass &parent, const std::string &name) : NXDataSet(parent, name), m_n(0) {}
  /** Returns a pointer to the internal data buffer.
   *  @throw runtime_error exception if the data have not been loaded /
   * initialized.
   *  @return a pointer to the array of items
   */
  const T *operator()() const {
    if (m_data.empty())
      throw std::runtime_error("Attempt to read uninitialized data from " + path());
    return m_data.data();
  }

  T *operator()() {
    if (m_data.empty())
      throw std::runtime_error("Attempt to read uninitialized data from " + path());
    return m_data.data();
  }

  /** Returns the i-th value in the internal buffer
   *  @param i :: The linear index of the data element
   *  @throw runtime_error if the data have not been loaded / initialized.
   *  @throw range_error if the index is greater than the buffer size.
   *  @return A reference to the value
   */
  const T &operator[](int i) const {
    if (m_data.empty())
      throw std::runtime_error("Attempt to read uninitialized data from " + path());
    if (i < 0 || i >= m_n)
      rangeError();
    return m_data[i];
  }

  T &operator[](int i) { return const_cast<T &>(static_cast<const NXDataSetTyped &>(*this)[i]); }
  /** Returns a value assuming the data is a two-dimensional array
   *  @param i :: The index along dim0()
   *  @param j :: The index along dim1()
   *  @throw runtime_error if the data have not been loaded / initialized.
   *  @throw range_error if the indeces point outside the buffer.
   *  @return A reference to the value
   */
  const T &operator()(int i, int j) const { return this->operator[](i * dim1() + j); }
  T &operator()(int i, int j) { return const_cast<T &>(static_cast<const NXDataSetTyped &>(*this)(i, j)); }
  /** Returns a value assuming the data is a tree-dimensional array
   *  @param i :: The index along dim0()
   *  @param j :: The index along dim1()
   *  @param k :: The index along dim2()
   *  @throw runtime_error if the data have not been loaded / initialized.
   *  @throw range_error if the indeces point outside the buffer.
   *  @return A reference to the value
   */
  const T &operator()(int i, int j, int k) const { return this->operator[]((i * dim1() + j) * dim2() + k); }
  T &operator()(int i, int j, int k) { return const_cast<T &>(static_cast<const NXDataSetTyped &>(*this)(i, j, k)); }

  /// Returns a the internal buffer
  container_T<T> &vecBuffer() { return m_data; }
  /// Returns the size of the data buffer
  int size() const { return m_n; }
  /**  Implementation of the virtual NXDataSet::load(...) method. Internally the
   * data are stored as a 1d array. If the data are loaded in chunks the newly read in
   * data replace the old ones. The actual rank of the loaded
   * data is equal or less than the rank of the dataset (returned by rank()
   * method).
   * Reads in all data.
   */
  void load() override {
    if (rank() > 4) {
      throw std::runtime_error("Cannot load dataset of rank greater than 4");
    }
    if (rank() == 4) {
      int n = dim0() * dim1() * dim2() * dim3();
      alloc(n);
      getData(m_data.data());
    } else if (rank() == 3) {
      int n = dim0() * dim1() * dim2();
      alloc(n);
      getData(m_data.data());
    } else if (rank() == 2) {
      int n = dim0() * dim1();
      alloc(n);
      getData(m_data.data());
    } else if (rank() == 1) {
      int n = dim0();
      alloc(n);
      getData(m_data.data());
    }
  }

private:
  /** Allocates memory for the data buffer
   *  @param n :: The number of elements to allocate.
   */
  void alloc(int n) {
    if (n <= 0) {
      throw std::runtime_error("Attempt to load from an empty dataset " + path());
    }
    try {
      if (m_n != n) {
        m_data.resize(n);
        m_n = n;
      }
    } catch (...) {
      std::stringstream ostr;
      ostr << "Cannot allocate " << (n * sizeof(T)) << " bytes of memory to load the data";
      throw std::runtime_error(ostr.str());
    }
  }
  /// A shortcut to "throw std::range_error("Nexus dataset range error");"
  void rangeError() const { throw std::range_error("Nexus dataset range error"); }
  // We cannot use an STL vector due to the dreaded std::vector<bool>
  container_T<T> m_data; ///< The data buffer
  int m_size[4];         ///< The sizes of the loaded data
  int m_n;               ///< The buffer size
};

/// The integer dataset type
using NXInt = NXDataSetTyped<int>;
/// The float dataset type
using NXFloat = NXDataSetTyped<float>;
/// The double dataset type
using NXDouble = NXDataSetTyped<double>;
/// The char dataset type
using NXChar = NXDataSetTyped<char>;
/// The size_t dataset type
using NXSize = NXDataSetTyped<std::size_t>;
/// The size_t dataset type
using NXUInt = NXDataSetTyped<unsigned int>;

//-------------------- classes --------------------------//

/**  The base class for a Nexus class (group). A Nexus class can contain
 * datasets and other Nexus classes.
 *   The NeXus file format (www.nexusformat.org) specifies the content of the
 * Nexus classes.
 *   Derived classes have specialized methods for creating classes and datasets
 * specific for the particular Nexus class.
 *   NXClass is a conctrete C++ class so arbitrary, non-standard Nexus classes
 * (groups) can be created and loaded from
 *   NeXus files.
 */
class MANTID_LEGACYNEXUS_DLL NXClass : public NXObject {
  friend class NXRoot;

public:
  /**  Constructor.
   *   @param parent :: The parent Nexus class. In terms of HDF it is the group
   * containing the NXClass.
   *   @param name :: The name of the NXClass relative to its parent
   */
  NXClass(const NXClass &parent, const std::string &name);
  /// The NX class identifier
  std::string NX_class() const override { return "NXClass"; }
  /**  Returns the class information about the next entry (class or dataset) in
   * this class.
   */
  NXClassInfo getNextEntry();
  /// Creates a new object in the NeXus file at path path.
  // virtual void make(const std::string& path) = 0;
  /// Resets the current position for getNextEntry() to the beginning
  void reset();
  /**  Templated method for creating derived NX classes. It also opens the
   * created class.
   *   @param name :: The name of the class
   *   @tparam NX Concrete Nexus class
   *   @return The new object
   */
  template <class NX> NX openNXClass(const std::string &name) const {
    NX nxc(*this, name);
    nxc.open();
    return nxc;
  }

  /**  Creates and opens an arbitrary (non-standard) class (group).
   *   @param name :: The name of the class.
   *   @return The opened NXClass
   */
  NXClass openNXGroup(const std::string &name) const { return openNXClass<NXClass>(name); }

  /**  Templated method for creating datasets. It also opens the created set.
   *   @param name :: The name of the dataset
   *   @tparam T The type of the data (int, double, ...).
   *   @return The new object
   */
  template <class T> NXDataSetTyped<T> openNXDataSet(const std::string &name) const {
    NXDataSetTyped<T> data(*this, name);
    data.open();
    return data;
  }

  /**  Creates and opens an integer dataset
   *   @param name :: The name of the dataset
   *   @return The int
   */
  NXInt openNXInt(const std::string &name) const { return openNXDataSet<int>(name); }
  /**  Creates and opens a float dataset
   *   @param name :: The name of the dataset
   *   @return The float
   */
  NXFloat openNXFloat(const std::string &name) const { return openNXDataSet<float>(name); }
  /**  Creates and opens a double dataset
   *   @param name :: The name of the dataset
   *   @return The double
   */
  NXDouble openNXDouble(const std::string &name) const { return openNXDataSet<double>(name); }
  /**  Creates and opens a char dataset
   *   @param name :: The name of the dataset
   *   @return The char
   */
  NXChar openNXChar(const std::string &name) const { return openNXDataSet<char>(name); }
  /**  Creates and opens a size_t dataset
   *   @param name :: The name of the dataset
   *   @return The size_t
   */
  NXSize openNXSize(const std::string &name) const { return openNXDataSet<std::size_t>(name); }
  /**  Returns a string
   *   @param name :: The name of the NXChar dataset
   *   @return The string
   */
  std::string getString(const std::string &name) const;
  /**  Returns a double
   *   @param name :: The name of the NXDouble dataset
   *   @return The double
   */
  double getDouble(const std::string &name) const;
  /**  Returns a float
   *   @param name :: The name of the NXFloat dataset
   *   @return The float
   */
  float getFloat(const std::string &name) const;
  /**  Returns a int
   *   @param name :: The name of the NXInt dataset
   *   @return The int
   */
  int getInt(const std::string &name) const;

  /// Returns a list of all classes (or groups) in this NXClass
  std::vector<NXClassInfo> &groups() const { return *m_groups; }
  /// Returns whether an individual group (or group) is present
  bool containsGroup(const std::string &query) const;
  /// Returns a list of all datasets in this NXClass
  std::vector<NXInfo> &datasets() const { return *m_datasets; }
  /** Returns NXInfo for a dataset
   *  @param name :: The name of the dataset
   *  @return NXInfo::stat is set to NXstatus::NX_ERROR if the dataset does not exist
   */
  NXInfo getDataSetInfo(const std::string &name) const;
  /// Returns whether an individual dataset is present
  bool containsDataSet(const std::string &query) const;
  /// Opens this NXClass using NXopengrouppath. Can be slow (or is slow)
  void open();

protected:
  std::shared_ptr<std::vector<NXClassInfo>> m_groups; ///< Holds info about the child NXClasses
  std::shared_ptr<std::vector<NXInfo>> m_datasets;    ///< Holds info about the datasets in this NXClass
  void readAllInfo();                                 ///< Fills in m_groups and m_datasets.
  void clear();                                       ///< Deletes content of m_groups and m_datasets
private:
  /// Pricate constructor.
  NXClass() : NXObject() { clear(); }
};

//-------------------- main classes -------------------------------//

/**  Implements NXdata Nexus class.
 */
class MANTID_LEGACYNEXUS_DLL NXData : public NXClass {
public:
  /**  Constructor.
   *   @param parent :: The parent Nexus class. In terms of HDF it is the group
   * containing the NXClass.
   *   @param name :: The name of the NXClass relative to its parent
   */
  NXData(const NXClass &parent, const std::string &name);
  /// Nexus class id
  std::string NX_class() const override { return "NXdata"; }
  /**  Opens the dataset within this NXData with signal=1 attribute.
   */
  template <typename T> NXDataSetTyped<T> openData() {
    for (std::vector<NXInfo>::const_iterator it = datasets().begin(); it != datasets().end(); ++it) {
      NXDataSet dset(*this, it->nxname);
      dset.open();
      // std::cerr << "NXData signal of " << it->nxname << " = " <<
      // dset.attributes("signal") << "\n";
      if (dset.attributes("signal") == "1") {
        return openNXDataSet<T>(it->nxname);
      }
    }
    // You failed to find the signal.
    // So try to just open the "data" entry directly
    return openNXDataSet<T>("data");
    // throw std::runtime_error("NXData does not seem to contain the data");
    // return NXDataSetTyped<T>(*this,"");
  }
  /// Opens data of double type
  NXDouble openDoubleData() { return openData<double>(); }
  /// Opens data of float type
  NXFloat openFloatData() { return openData<float>(); }
  /// Opens data of int type
  NXInt openIntData() { return openData<int>(); }
  /// Opens data of size type
  NXSize openSizeData() { return openData<std::size_t>(); }
  /// Opens data of unsigned int type
  NXUInt openUIntData() { return openData<unsigned int>(); }
};

/**  Implements NXdetector Nexus class.
 */
class MANTID_LEGACYNEXUS_DLL NXDetector : public NXClass {
public:
  /**  Constructor.
   *   @param parent :: The parent Nexus class. In terms of HDF it is the group
   * containing the NXClass.
   *   @param name :: The name of the NXClass relative to its parent
   */
  NXDetector(const NXClass &parent, const std::string &name) : NXClass(parent, name) {}
  /// Nexus class id
  std::string NX_class() const override { return "NXdetector"; }
  /// Opens the dataset containing pixel distances
  NXFloat openDistance() { return openNXFloat("distance"); }
  /// Opens the dataset containing pixel azimuthal angles
  NXFloat openAzimuthalAngle() { return openNXFloat("azimuthal_angle"); }
  /// Opens the dataset containing pixel polar angles
  NXFloat openPolarAngle() { return openNXFloat("polar_angle"); }
};

/**  Implements NXinstrument Nexus class.
 */
class MANTID_LEGACYNEXUS_DLL NXInstrument : public NXClass {
public:
  /**  Constructor.
   *   @param parent :: The parent Nexus class. In terms of HDF it is the group
   * containing the NXClass.
   *   @param name :: The name of the NXClass relative to its parent
   */
  NXInstrument(const NXClass &parent, const std::string &name) : NXClass(parent, name) {}
  /// Nexus class id
  std::string NX_class() const override { return "NXinstrument"; }
  /**  Opens a NXDetector
   *   @param name :: The name of the class
   *   @return The detector
   */
  NXDetector openNXDetector(const std::string &name) { return openNXClass<NXDetector>(name); }
};

/**  Implements NXentry Nexus class.
 */
class MANTID_LEGACYNEXUS_DLL NXEntry : public NXClass {
public:
  /**  Constructor.
   *   @param parent :: The parent Nexus class. In terms of HDF it is the group
   * containing the NXClass.
   *   @param name :: The name of the NXClass relative to its parent
   */
  NXEntry(const NXClass &parent, const std::string &name) : NXClass(parent, name) {}
  /// Nexus class id
  std::string NX_class() const override { return "NXentry"; }
  /**  Opens a NXData
   *   @param name :: The name of the class
   *   @return the nxdata entry
   */
  NXData openNXData(const std::string &name) const { return openNXClass<NXData>(name); }
  /**  Opens a NXInstrument
   *   @param name :: The name of the class
   *   @return the instrument
   */
  NXInstrument openNXInstrument(const std::string &name) const { return openNXClass<NXInstrument>(name); }
};

/**  Implements NXroot Nexus class.
 */
class MANTID_LEGACYNEXUS_DLL NXRoot : public NXClass {
public:
  // Constructor
  NXRoot(std::string fname);
  /// Destructor
  ~NXRoot() override;
  /// Return the NX class for a class (HDF group) or "SDS" for a data set;
  std::string NX_class() const override { return "NXroot"; }
  /**  Opens an entry -- a topmost Nexus class
   *   @param name :: The name of the entry
   *   @return the entry
   */
  NXEntry openEntry(const std::string &name) { return openNXClass<NXEntry>(name); }
  NXEntry openFirstEntry();

private:
  const std::string m_filename; ///< The file name
};

} // namespace LegacyNexus
} // namespace Mantid
