// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
#include <cstdint>
#include <vector>
#include <xdiag/bits/bitvector.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::symmetries {

template <typename enumeration_tt> class RepresentativeTable {
public:
  using enumeration_t = enumeration_tt;
  using bit_t = typename enumeration_t::bit_t;

  RepresentativeTable() = default;
  RepresentativeTable(enumeration_t const &enumeration,
                      SitePermutation const &action,
                      Representation const &irrep);

  inline int64_t representative(int64_t idx) const {
    return representative_[representative_index_[idx]];
  }
  inline int64_t index_of_representative(int64_t idx) const {
    return representative_index_[idx];
  }
  inline int64_t symmetry(int64_t idx) const {
    return representative_symmetry_[idx];
  }
  inline int64_t norm(int64_t idx) const {
    return norm_[representative_norm_index_[idx]];
  }
  int64_t size() const;

  bool operator==(RepresentativeTable<enumeration_t> const &rhs) const;
  bool operator!=(RepresentativeTable<enumeration_t> const &rhs) const;

private:
  bits::BitVector<bit_t> representative_;
  bits::BitVector<bit_t> representative_index_;
  bits::BitVector<bit_t> representative_symmetry_;
  bits::BitVector<bit_t> representative_norm_index_;
  std::vector<double> norm_;
};

} // namespace xdiag::symmetries
