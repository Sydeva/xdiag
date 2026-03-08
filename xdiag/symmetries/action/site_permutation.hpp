// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class SitePermutation {
public:
  XDIAG_API SitePermutation() = default;
  XDIAG_API explicit SitePermutation(PermutationGroup const &group);

  XDIAG_API int64_t size() const;
  XDIAG_API int64_t nsites() const;
  template <typename bit_t> XDIAG_API bit_t apply(int64_t idx, bit_t bits);

  XDIAG_API bool operator==(SitePermutation const &rhs) const;
  XDIAG_API bool operator!=(SitePermutation const &rhs) const;

private:
  PermutationGroup group_;
};

} // namespace xdiag
