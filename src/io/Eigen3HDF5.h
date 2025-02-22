/**
 * @file Eigen3HDF5.h
 *
 * @date Jun 8, 2016
 * @author Jan Unsleber
 *
 * Taken and modified from:
 * https://github.com/garrison/eigen3-hdf5/commits/master
 * Tag: 2c782414251e75a2de9b0441c349f5f18fe929a2
 *
 * @copyright
 *  SPECIAL LICENSE FOR THIS FILE ONLY
 *
 *  The MIT License (MIT)
 *
 *  Copyright (c) 2013 James R. Garrison
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

#ifndef _EIGEN3HDF5_H_
#define _EIGEN3HDF5_H_
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <H5Cpp.h>
#include <hdf5.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace Serenity {
/**
 * @namespace EigenHDF5 EigenHDF5.h
 * @brief Namespace for all Eigen3 to HDF5 bindings.
 */
namespace EigenHDF5 {

inline H5::DataSpace getSpace(const H5::DataSet& ds) {
  hid_t id2 = ds.getId();
  hid_t myspace = H5Dget_space(id2);
  H5::DataSpace origspace(myspace);
  H5Sclose(myspace);
  return origspace;
}

template<typename T>
struct DatatypeSpecialization;

// floating-point types

template<>
struct DatatypeSpecialization<float> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_FLOAT;
  }
};

template<>
struct DatatypeSpecialization<double> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_DOUBLE;
  }
};

template<>
struct DatatypeSpecialization<long double> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_LDOUBLE;
  }
};

// integer types

template<>
struct DatatypeSpecialization<short> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_SHORT;
  }
};

template<>
struct DatatypeSpecialization<unsigned short> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_USHORT;
  }
};

template<>
struct DatatypeSpecialization<int> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_INT;
  }
};

template<>
struct DatatypeSpecialization<unsigned int> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_UINT;
  }
};

template<>
struct DatatypeSpecialization<long> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_LONG;
  }
};

template<>
struct DatatypeSpecialization<unsigned long> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_ULONG;
  }
};

template<>
struct DatatypeSpecialization<long long> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_LLONG;
  }
};

template<>
struct DatatypeSpecialization<unsigned long long> {
  static inline const H5::DataType* get(void) {
    return &H5::PredType::NATIVE_ULLONG;
  }
};

// complex types
//
// inspired by http://www.mail-archive.com/hdf-forum@hdfgroup.org/msg00759.html

template<typename T>
class ComplexH5Type : public H5::CompType {
 public:
  ComplexH5Type(void) : CompType(sizeof(std::complex<T>)) {
    const H5::DataType* const datatype = DatatypeSpecialization<T>::get();
    assert(datatype->getSize() == sizeof(T));
    // If we call the members "r" and "i", h5py interprets the
    // structure correctly as complex numbers.
    this->insertMember(std::string("r"), 0, *datatype);
    this->insertMember(std::string("i"), sizeof(T), *datatype);
    this->pack();
  }

  static const ComplexH5Type<T>* get_singleton(void) {
    // NOTE: constructing this could be a race condition
    static ComplexH5Type<T> singleton;
    return &singleton;
  }
};

template<typename T>
struct DatatypeSpecialization<std::complex<T>> {
  static inline const H5::DataType* get(void) {
    return ComplexH5Type<T>::get_singleton();
  }
};

// string types, to be used mainly for attributes

template<>
struct DatatypeSpecialization<const char*> {
  static inline const H5::DataType* get(void) {
    static const H5::StrType strtype(0, H5T_VARIABLE);
    return &strtype;
  }
};

template<>
struct DatatypeSpecialization<char*> {
  static inline const H5::DataType* get(void) {
    static const H5::StrType strtype(0, H5T_VARIABLE);
    return &strtype;
  }
};
template<>
struct DatatypeSpecialization<std::string> {
  static inline const H5::DataType* get(void) {
    static const H5::StrType strtype(0, H5T_VARIABLE);
    return &strtype;
  }
};

// XXX: for some unknown reason the following two functions segfault if
// H5T_VARIABLE is used.  The passed strings should still be null-terminated,
// so this is a bit worrisome.

template<std::size_t N>
struct DatatypeSpecialization<const char[N]> {
  static inline const H5::DataType* get(void) {
    static const H5::StrType strtype(0, N);
    return &strtype;
  }
};

template<std::size_t N>
struct DatatypeSpecialization<char[N]> {
  static inline const H5::DataType* get(void) {
    static const H5::StrType strtype(0, N);
    return &strtype;
  }
};

namespace internal {
template<typename Derived>
H5::DataSpace create_dataspace(const Eigen::EigenBase<Derived>& mat) {
  const std::size_t dimensions_size = 2;
  const hsize_t dimensions[dimensions_size] = {static_cast<hsize_t>(mat.rows()), static_cast<hsize_t>(mat.cols())};
  return H5::DataSpace(dimensions_size, dimensions);
}

template<typename Derived>
bool write_rowmat(const Eigen::EigenBase<Derived>& mat, const H5::DataType* const datatype, H5::DataSet* dataset,
                  const H5::DataSpace* dspace) {
  if (mat.derived().innerStride() != 1) {
    // inner stride != 1 is an edge case this function does not (yet) handle. (I think it
    // could by using the inner stride as the first element of mstride below. But I do
    // not have a test case to try it out, so just return false for now.)
    return false;
  }

  assert(mat.rows() >= 0);
  assert(mat.cols() >= 0);
  assert(mat.derived().outerStride() >= 0);
  hsize_t rows = hsize_t(mat.rows());
  hsize_t cols = hsize_t(mat.cols());
  hsize_t stride = hsize_t(mat.derived().outerStride());

  // slab params for the file data
  hsize_t fstride[2] = {1, cols};

  // slab params for the memory data
  hsize_t mstride[2] = {1, stride};

  // slab params for both file and memory data
  hsize_t count[2] = {1, 1};
  hsize_t block[2] = {rows, cols};
  hsize_t start[2] = {0, 0};

  // memory dataspace
  hsize_t mdim[2] = {rows, stride};
  H5::DataSpace mspace(2, mdim);

  dspace->selectHyperslab(H5S_SELECT_SET, count, start, fstride, block);
  mspace.selectHyperslab(H5S_SELECT_SET, count, start, mstride, block);
  dataset->write(mat.derived().data(), *datatype, mspace, *dspace);

  return true;
}

template<typename Derived>
bool write_colmat(const Eigen::EigenBase<Derived>& mat, const H5::DataType* const datatype, H5::DataSet* dataset,
                  const H5::DataSpace* dspace) {
  if (mat.derived().innerStride() != 1) {
    // inner stride != 1 is an edge case this function does not (yet) handle. (I think it
    // could by using the inner stride as the first element of mstride below. But I do
    // not have a test case to try it out, so just return false for now.)
    return false;
  }

  assert(mat.rows() >= 0);
  assert(mat.cols() >= 0);
  assert(mat.derived().outerStride() >= 0);
  hsize_t rows = hsize_t(mat.rows());
  hsize_t cols = hsize_t(mat.cols());
  hsize_t stride = hsize_t(mat.derived().outerStride());

  // slab params for the file data
  hsize_t fstride[2] = {1, cols};
  hsize_t fcount[2] = {1, 1};
  hsize_t fblock[2] = {1, cols};

  // slab params for the memory data
  hsize_t mstride[2] = {stride, 1};
  hsize_t mcount[2] = {1, 1};
  hsize_t mblock[2] = {cols, 1};

  // memory dataspace
  hsize_t mdim[2] = {cols, stride};
  H5::DataSpace mspace(2, mdim);

  // transpose the column major data in memory to the row major data in the file by
  // writing one row slab at a time.
  for (hsize_t i = 0; i < rows; i++) {
    hsize_t fstart[2] = {i, 0};
    hsize_t mstart[2] = {0, i};
    dspace->selectHyperslab(H5S_SELECT_SET, fcount, fstart, fstride, fblock);
    mspace.selectHyperslab(H5S_SELECT_SET, mcount, mstart, mstride, mblock);
    dataset->write(mat.derived().data(), *datatype, mspace, *dspace);
  }
  return true;
}
} // namespace internal

template<typename T>
void save_scalar_attribute(const H5::H5Object& h5obj, const std::string& name, const T& value) {
  const H5::DataType* const datatype = DatatypeSpecialization<T>::get();
  H5::DataSpace dataspace(H5S_SCALAR);
  H5::Attribute att = h5obj.createAttribute(name.c_str(), *datatype, dataspace);
  att.write(*datatype, &value);
}

template<>
inline void save_scalar_attribute(const H5::H5Object& h5obj, const std::string& name, const std::string& value) {
  save_scalar_attribute(h5obj, name, value.c_str());
}

/**
 * @brief Loads an attribute from file. The attribute is currently
 *        only saved as string and gets converted into typename T
 *        via stringstream.
 *        TODO Find a possibility to read other types directly.
 * @param h5obj The HDF5 file containing the attribute.
 * @param name The name of the attribute within the file
 * @param value The value to be filled with the attributes value.
 */
template<typename T>
void load_scalar_attribute(const H5::H5Object& h5obj, const std::string& name, T& value) {
  H5::Attribute attr = h5obj.openAttribute(name);
  DatatypeSpecialization<char*> mem_type;
  std::string buf;
  attr.read(*(mem_type.get()), buf);
  std::stringstream linestream(buf);
  linestream >> value;
}

template<typename Derived>
void replace(H5::H5Location& h5group, const std::string& name, const Eigen::EigenBase<Derived>& mat) {
  typedef typename Derived::Scalar Scalar;
  const H5::DataType* const datatype = DatatypeSpecialization<Scalar>::get();
  H5::DataSet dataset = h5group.openDataSet(name.c_str());
  const H5::DataSpace dataspace = dataset.getSpace();

  bool written = false; // flag will be true when the data has been written
  if (mat.derived().Flags & Eigen::RowMajor) {
    written = internal::write_rowmat(mat, datatype, &dataset, &dataspace);
  }
  else {
    written = internal::write_colmat(mat, datatype, &dataset, &dataspace);
  }

  if (!written) {
    // data has not yet been written, so there is nothing else to try but copy the input
    // matrix to a row major matrix and write it.
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row_major_mat(mat);
    dataset.write(row_major_mat.data(), *datatype);
  }
}

// see http://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html

template<typename Derived>
void save(H5::H5Location& h5group, const std::string& name, const Eigen::EigenBase<Derived>& mat,
          const H5::DSetCreatPropList& plist = H5::DSetCreatPropList::DEFAULT) {
  typedef typename Derived::Scalar Scalar;
  const H5::DataType* const datatype = DatatypeSpecialization<Scalar>::get();
  const H5::DataSpace dataspace = internal::create_dataspace(mat);
  H5::DataSet dataset = h5group.createDataSet(name.c_str(), *datatype, dataspace, plist);

  bool written = false; // flag will be true when the data has been written
  if (mat.derived().Flags & Eigen::RowMajor) {
    written = internal::write_rowmat(mat, datatype, &dataset, &dataspace);
  }
  else {
    written = internal::write_colmat(mat, datatype, &dataset, &dataspace);
  }

  if (!written) {
    // data has not yet been written, so there is nothing else to try but copy the input
    // matrix to a row major matrix and write it.
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row_major_mat(mat);
    dataset.write(row_major_mat.data(), *datatype);
  }
}

template<typename Derived>
void save_attribute(const H5::H5Object& h5obj, const std::string& name, const Eigen::EigenBase<Derived>& mat) {
  typedef typename Derived::Scalar Scalar;
  const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row_major_mat(mat);
  const H5::DataSpace dataspace = internal::create_dataspace(mat);
  const H5::DataType* const datatype = DatatypeSpecialization<Scalar>::get();
  H5::Attribute dataset = h5obj.createAttribute(name, *datatype, dataspace);
  dataset.write(*datatype, row_major_mat.data());
}

namespace internal {
// H5::Attribute and H5::DataSet both have similar API's, and although they
// share a common base class, the relevant methods are not virtual.  Worst
// of all, they take their arguments in different orders!

template<typename Scalar>
inline void read_data(const H5::DataSet& dataset, Scalar* data, const H5::DataType& datatype) {
  dataset.read(data, datatype);
}

template<typename Scalar>
inline void read_data(const H5::Attribute& dataset, Scalar* data, const H5::DataType& datatype) {
  dataset.read(datatype, data);
}

// read a column major attribute; I do not know if there is an hdf routine to read an
// attribute hyperslab, so I take the lazy way out: just read the conventional hdf
// row major data and let eigen copy it into mat.
template<typename Derived>
bool read_colmat(const Eigen::DenseBase<Derived>& mat, const H5::DataType* const datatype, const H5::Attribute& dataset) {
  typename Derived::Index rows = mat.rows();
  typename Derived::Index cols = mat.cols();
  typename Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> temp(rows, cols);
  internal::read_data(dataset, temp.data(), *datatype);
  const_cast<Eigen::DenseBase<Derived>&>(mat) = temp;
  return true;
}

template<typename Derived>
bool read_colmat(const Eigen::DenseBase<Derived>& mat, const H5::DataType* const datatype, const H5::DataSet& dataset) {
  if (mat.derived().innerStride() != 1) {
    // inner stride != 1 is an edge case this function does not (yet) handle. (I think it
    // could by using the inner stride as the first element of mstride below. But I do
    // not have a test case to try it out, so just return false for now.)
    return false;
  }

  assert(mat.rows() >= 0);
  assert(mat.cols() >= 0);
  assert(mat.derived().outerStride() >= 0);
  hsize_t rows = hsize_t(mat.rows());
  hsize_t cols = hsize_t(mat.cols());
  hsize_t stride = hsize_t(mat.derived().outerStride());

  if (stride != rows) {
    // this function does not (yet) read into a mat that has a different stride than the
    // dataset.
    return false;
  }

  // slab params for the file data
  hsize_t fstride[2] = {1, cols};
  hsize_t fcount[2] = {1, 1};
  hsize_t fblock[2] = {1, cols};

  // file dataspace
  hsize_t fdim[2] = {rows, cols};
  H5::DataSpace fspace(2, fdim);

  // slab params for the memory data
  hsize_t mstride[2] = {stride, 1};
  hsize_t mcount[2] = {1, 1};
  hsize_t mblock[2] = {cols, 1};

  // memory dataspace
  hsize_t mdim[2] = {cols, stride};
  H5::DataSpace mspace(2, mdim);

  // transpose the column major data in memory to the row major data in the file by
  // writing one row slab at a time.
  for (hsize_t i = 0; i < rows; i++) {
    hsize_t fstart[2] = {i, 0};
    hsize_t mstart[2] = {0, i};
    fspace.selectHyperslab(H5S_SELECT_SET, fcount, fstart, fstride, fblock);
    mspace.selectHyperslab(H5S_SELECT_SET, mcount, mstart, mstride, mblock);
    dataset.read(const_cast<Eigen::DenseBase<Derived>&>(mat).derived().data(), *datatype, mspace, fspace);
  }
  return true;
}

template<typename Derived, typename DataSet>
void _load(const DataSet& dataset, const Eigen::DenseBase<Derived>& mat) {
  typedef typename Derived::Scalar Scalar;
  const H5::DataSpace dataspace = dataset.getSpace();
  const std::size_t ndims = dataspace.getSimpleExtentNdims();
  assert(ndims > 0);
  const std::size_t dimensions_size = 2;
  hsize_t dimensions[dimensions_size];
  dimensions[1] = 1; // in case it's 1D
  if (ndims > dimensions_size) {
    throw SerenityError("HDF5 array has too many dimensions.");
  }
  dataspace.getSimpleExtentDims(dimensions);
  const hsize_t rows = dimensions[0], cols = dimensions[1];
  const H5::DataType* const datatype = DatatypeSpecialization<Scalar>::get();
  Eigen::DenseBase<Derived>& mat_ = const_cast<Eigen::DenseBase<Derived>&>(mat);
  mat_.derived().resize(rows, cols);
  bool written = false;
  bool isRowMajor = mat.Flags & Eigen::RowMajor;
  if (isRowMajor || dimensions[0] == 1 || dimensions[1] == 1) {
    // mat is already row major
    typename Derived::Index istride = mat_.derived().outerStride();
    assert(istride >= 0);
    hsize_t stride = istride >= 0 ? istride : 0;
    if (stride == cols || (stride == rows && cols == 1)) {
      // mat has natural stride, so read directly into its data block
      read_data(dataset, mat_.derived().data(), *datatype);
      written = true;
    }
  }
  else {
    // colmajor flag is 0 so the assert needs to check that mat is not rowmajor.
    assert(!(mat.Flags & Eigen::RowMajor));

    written = read_colmat(mat_, datatype, dataset);
  }

  if (!written) {
    // dataset has not been loaded directly into mat_, so as a last resort read it into a
    // temp and copy it to mat_. (Should only need to do this when the mat_ to be loaded
    // into has an unnatural stride.)
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> temp(rows, cols);
    internal::read_data(dataset, temp.data(), *datatype);
    mat_ = temp;
    written = true;
  }
}
} // namespace internal

template<typename Derived>
void load(const H5::H5Location& h5group, const std::string& name, const Eigen::DenseBase<Derived>& mat) {
  const H5::DataSet dataset = h5group.openDataSet(name.c_str());
  internal::_load(dataset, mat);
}

template<typename Derived>
void load_attribute(const H5::H5Object& h5obj, const std::string& name, const Eigen::DenseBase<Derived>& mat) {
  const H5::Attribute dataset = h5obj.openAttribute(name.c_str());
  internal::_load(dataset, mat);
}

// Caution! Only use this function if you are sure the target and source object have the
// same type (double based) and dimensions
//
// This function is used to drastically reduce overhead for I/O bottlenecked applications
template<typename Derived>
void loadDirectly(const H5::DataSet& dataset, Eigen::DenseBase<Derived>& mat) {
  const H5::DataType* const doubletype = DatatypeSpecialization<double>::get();
  internal::read_data(dataset, mat.derived().data(), *doubletype);
}

template<typename Scalar>
class MyTriplet : public Eigen::Triplet<Scalar> {
 public:
  MyTriplet(void) : Eigen::Triplet<Scalar>() {
  }

  MyTriplet(const unsigned int& i, const unsigned int& j, const Scalar& v = Scalar(0))
    : Eigen::Triplet<Scalar>(i, j, v) {
  }

  static std::size_t offsetof_row(void) {
    return offsetof(MyTriplet<Scalar>, m_row);
  }
  static std::size_t offsetof_col(void) {
    return offsetof(MyTriplet<Scalar>, m_col);
  }
  static std::size_t offsetof_value(void) {
    return offsetof(MyTriplet<Scalar>, m_value);
  }
};

template<typename T>
class SparseH5Type : public H5::CompType {
 public:
  SparseH5Type(void) : CompType(sizeof(MyTriplet<T>)) {
    const H5::DataType* const datatypei = DatatypeSpecialization<unsigned int>::get();
    const H5::DataType* const datatype = DatatypeSpecialization<T>::get();
    assert(datatype->getSize() == sizeof(T));
    this->insertMember(std::string("r"), MyTriplet<T>::offsetof_row(), *datatypei);
    this->insertMember(std::string("c"), MyTriplet<T>::offsetof_col(), *datatypei);
    this->insertMember(std::string("v"), MyTriplet<T>::offsetof_value(), *datatype);
    this->pack();
  }

  static const SparseH5Type<T>* get_singleton(void) {
    // NOTE: constructing this could be a race condition
    static SparseH5Type<T> singleton;
    return &singleton;
  }
};

template<typename SparseMatrixType>
void save_sparse(H5::H5Location& h5group, const std::string& name, const SparseMatrixType& mat,
                 const H5::DSetCreatPropList& plist = H5::DSetCreatPropList::DEFAULT) {
  typedef typename SparseMatrixType::Scalar Scalar;
  // save the actual sparse matrix
  std::vector<MyTriplet<Scalar>> data;
  data.reserve(mat.nonZeros());
  for (int k = 0; k < mat.outerSize(); ++k) {
    for (typename SparseMatrixType::InnerIterator it(mat, k); it; ++it) {
      if (it.value() != Scalar(0))
        data.push_back(MyTriplet<Scalar>(it.row(), it.col(), it.value()));
    }
  }
  const hsize_t nnz = data.size();
  const H5::DataSpace dataspace(1, &nnz);
  const H5::DataType* const datatype = SparseH5Type<Scalar>::get_singleton();
  H5::DataSet dataset = h5group.createDataSet(name, *datatype, dataspace, plist);
  dataset.write(data.data(), *datatype);
  // save the matrix's shape as an attribute
  Eigen::Matrix<typename SparseMatrixType::Index, 2, 1> shape;
  shape(0) = mat.rows();
  shape(1) = mat.cols();
  save_attribute(dataset, "shape", shape);
}

template<typename SparseMatrixType>
void load_sparse(const H5::H5Location& h5group, const std::string& name, SparseMatrixType& mat) {
  typedef typename SparseMatrixType::Scalar Scalar;
  const H5::DataSet dataset = h5group.openDataSet(name);
  const H5::DataSpace dataspace = getSpace(dataset);
  const std::size_t ndims = dataspace.getSimpleExtentNdims();
  if (ndims != 1) {
    throw SerenityError("HDF5 array has incorrect number of dimensions to represent a sparse matrix.");
  }
  Eigen::Matrix<typename SparseMatrixType::Index, 2, 1> shape;
  load_attribute(dataset, "shape", shape);
  hsize_t nnz;
  dataspace.getSimpleExtentDims(&nnz); // assumes ndims == 1 in the data representation
  const H5::DataType* const datatype = SparseH5Type<Scalar>::get_singleton();
  std::vector<MyTriplet<Scalar>> data(nnz);
  dataset.read(data.data(), *datatype, dataspace);
  mat.resize(shape(0), shape(1)); // NOTE: this also clears all existing values
  mat.setFromTriplets(data.begin(), data.end());
}

} // namespace EigenHDF5
} // namespace Serenity
#endif
