/*
 *  Copyright (C) 2005 M.J. Zaki <zaki@cs.rpi.edu> Rensselaer Polytechnic Institute
 *  Written by parimi@cs.rpi.edu
 *  Updated by chaojv@cs.rpi.edu, alhasan@cs.rpi.edu, salems@cs.rpi.edu
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
 */
#ifndef _MSET_ISO_CHECK_H_
#define _MSET_ISO_CHECK_H_

#include "typedefs.h"

template<class PP, class MP, class ST, template<class, typename, typename, template <typename> class > class CC,
         template <typename> class ALLOC >
class pattern;

class no_edges;

template<class prop, class next_property>
class proplist;

/*
 * Method performs isomorphism checking for itemsets.
 * Also update the canonical code for the pattern.
 *
 * Method returns true to indicate that the pattern is canonical. For itemsets this is 
 * always true due to the manner in which we generate candidates. 
 */
template<typename PP, class MP, class PAT_ST, 
         template<class, typename, typename, template <typename> class > class CAN_CODE, 
         template <typename> class ALLOC >
bool check_isomorphism(pattern<MSET_PROP, MP, PAT_ST, CAN_CODE, ALLOC>* const& cand_pat) { 
  // NOTE:
  // The candidate generation process results in no isomorphism check for itemsets.

  // Update canonical code.
  cand_pat->_canonical_code.insert_vertex(cand_pat->rmost_vertex());

  cand_pat->_is_canonical = true;
  return true;
}

#endif
