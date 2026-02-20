// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "subsets.hpp"
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/ipow.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
Subsets<bit_t>::Subsets(int64_t n) try : n_(n), size_(utils::ipow(2, n)) {
  if (n < 0) {
    XDIAG_THROW("Error constructing Subsets: n<0");
  }
}
XDIAG_CATCH

template <class bit_t> int64_t Subsets<bit_t>::n() const { return n_; }
template <class bit_t> int64_t Subsets<bit_t>::size() const { return size_; }

template <class bit_t> SubsetsIterator<bit_t> Subsets<bit_t>::begin() const {
  return SubsetsIterator<bit_t>((int64_t)0);
}

template <class bit_t> SubsetsIterator<bit_t> Subsets<bit_t>::end() const {
  return SubsetsIterator<bit_t>((int64_t)size_);
}

template <class bit_t>
bool Subsets<bit_t>::operator==(Subsets<bit_t> const &rhs) const {
  return n_ == rhs.n_;
}

template <class bit_t>
bool Subsets<bit_t>::operator!=(Subsets<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class Subsets<uint16_t>;
template class Subsets<uint32_t>;
template class Subsets<uint64_t>;

template <class bit_t>
SubsetsIterator<bit_t>::SubsetsIterator(int64_t idx) : current_((bit_t)idx) {}

template <class bit_t>
bool SubsetsIterator<bit_t>::operator==(
    const SubsetsIterator<bit_t> &rhs) const {
  return current_ == rhs.current_;
}

template <class bit_t>
bool SubsetsIterator<bit_t>::operator!=(
    const SubsetsIterator<bit_t> &rhs) const {
  return !operator==(rhs);
}

template <class bit_t>
SubsetsIterator<bit_t> &SubsetsIterator<bit_t>::operator++() {
  ++current_;
  return *this;
}

template <class bit_t> bit_t SubsetsIterator<bit_t>::operator*() const {
  return current_;
}

template class SubsetsIterator<uint16_t>;
template class SubsetsIterator<uint32_t>;
template class SubsetsIterator<uint64_t>;

} // namespace xdiag::combinatorics
