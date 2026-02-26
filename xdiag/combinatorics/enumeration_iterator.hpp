// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once {

class EnumerationIterator {
public:
  virtual EnumerationIterator &operator++() = 0;
  virtual EnumerationIterator &operator*() const = 0;
  virtual bool operator==(EnumerationIterator const &rhs) const;
  virtual bool operator!=(EnumerationIterator const &rhs) const;
};

}
