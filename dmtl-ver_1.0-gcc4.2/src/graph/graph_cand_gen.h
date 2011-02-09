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
#ifndef _GRAPH_CAND_GEN_H
#define _GRAPH_CAND_GEN_H
#include "level_one_hmap.h"
#include "typedefs.h"
#include <algorithm>
#include <vector>
#include <string>

extern unsigned long int freq_pats_count;
extern bool print;

// candidate generation for graphs
template<typename PP, class MP, class PAT_ST, template<class, typename, typename, template <typename> class > class CC, template <typename> class ALLOC, class EDGE_MAP, class SM_TYPE>
void cand_gen(GRAPH_PATTERN* pat, const EDGE_MAP& emap, const int& minsup, pat_fam<GRAPH_PATTERN >& freq_pats, 
              count_support<GRAPH_PROP, V_Fk1_MINE_PROP, PAT_ST, CC, ALLOC, SM_TYPE>& cs) {

#ifdef PRINT
  cout<<"In call to cand_gen"<<endl;  
#endif
  typedef pat_fam<GRAPH_PATTERN> PAT_FAM;
 
  /// add back edges to get new candidates ///
  int back_signal = back_extensions(pat, emap, minsup, freq_pats, cs);

//	cout<<"back signal is "<<back_signal<<endl;
  /// add fwd edges to get new candidates ///
  int fwd_signal = fwd_extensions(pat, emap, minsup, freq_pats, cs); //signal value is one if it can be extended...0 o/w
//	cout<<"Forward extension done "<<endl;
//	cout<<"fwd signal is "<<fwd_signal<<endl;
	
	if((fwd_signal + back_signal) == 0) //Both the signals returned a zero which means that the pattern cannot be extended further...so it is max pat
	{
		cout<<"Reached a maximal pattern.."<<endl;
		cout<<pat<<endl;
	}



}//end cand_gen()


template<typename PP, class MP, class PAT_ST, template<class, typename, typename, template <typename> class> class CC, template <typename> class ALLOC, class EDGE_MAP, class SM_TYPE>
int back_extensions(GRAPH_PATTERN* pat, const EDGE_MAP& emap, const int& minsup, pat_fam<GRAPH_PATTERN >& freq_pats, 
                     count_support<GRAPH_PROP, V_Fk1_MINE_PROP, PAT_ST, CC, ALLOC, SM_TYPE>& cs) {

  if(pat->rmp_size()<3) // back extensions only for graphs with rmost path 
    return 0;             // with atleast two edges
#ifdef PRINT
  cout<<"Back extension for "<<pat<<endl;
#endif
  typedef pat_fam<GRAPH_PATTERN> PAT_FAM;

  const typename GRAPH_PATTERN::VERTEX_T& last_v=pat->rmost_vertex();
  const typename GRAPH_PATTERN::RMP_T& rmp=pat->rmost_path();
	pat->random_rmost_path(rmp);
//  const typename GRAPH_PATTERN::RMP_T& rmp=pat->random_rmost_path(); //Vertices on the rightmost path in a random order
  typename GRAPH_PATTERN::RMP_T::const_iterator rmp_it;
  typename GRAPH_PATTERN::EDGE_T e;
  GRAPH_PATTERN* edge=0;
  GRAPH_PATTERN* cand_pat=0;
  int rvid=pat->rmost_vid();
  int vid;

  // determine vid on rmp from where to start adding back edges
  rmp_it=rmp.end()-3;
  while(true) {
    if(rmp_it<rmp.begin())
      break;
    if(pat->get_out_edge(rvid, *rmp_it, e)) {
      rmp_it++;
      break;
    }
    rmp_it--;
  }

  if(rmp_it<rmp.begin())
    rmp_it++;

  while(rmp_it<rmp.end()-2) {
    vid=*rmp_it;
    const typename GRAPH_PATTERN::VERTEX_T& back_v=pat->label(vid);

    //For the directed case check if the edge exists first...
    if(!emap.has_label(last_v,back_v))
    {
	    rmp_it++;
	    continue;
    }
    const typename EDGE_MAP::LABELS& lbls=emap.get_labels(last_v, back_v);
    typename EDGE_MAP::CONST_LIT lit=lbls.begin();

    // get all possible labels for this back edge
    while(lit!=lbls.end()) {

      // construct new candidate with this label
      cand_pat=pat->clone();
      cand_pat->add_out_edge(rvid, vid, *lit);
      cand_pat->add_out_edge(vid, rvid, *lit);

      if(!check_isomorphism(cand_pat)) {
        delete cand_pat;
        cand_pat=0;
        lit++;
        continue;
      }

      // create the one edged-pattern
      edge=new GRAPH_PATTERN;
      if(last_v<back_v)
        make_edge(edge, last_v, back_v, *lit);
      else
        make_edge(edge, back_v, last_v, *lit);
#ifdef PRINT
      cout<<"Trying with back edge="<<edge<<endl;
      cout<<"New candidate="<<cand_pat<<endl;
#endif
      // count support
	//cout<<"Before count support"<<endl;
      cs.count(pat, edge, &cand_pat, minsup, 1);
	//cout<<"After count support"<<endl;

      delete edge;
      edge=0;

	if(!cand_pat)
	{
		lit++;
		continue;
	}
      // recursive call for frequent graph
      if(cand_pat->is_valid(minsup)) {
        //freq_pats.push_back(cand_pat);
	//cout<<"Recursive function call "<<endl;
        freq_pats_count++;
        if(print)
          cout << cand_pat;
        cand_gen(cand_pat, emap, minsup, freq_pats, cs);
        cs.delete_vat(cand_pat);
	return 1; // Candidate can be extended further
      }
      // else {
      delete cand_pat;
      cand_pat=0;
      // }

      lit++;
    }//end while lit

    rmp_it++;
  }//end while rmp_it
return 0; //Candidate cannot be extended further..
}//end back_extensions()


template<typename PP, class MP, class PAT_ST, template<typename, typename, typename, template <typename> class > class CC, template <typename> class ALLOC, class EDGE_MAP, class SM_TYPE>
int fwd_extensions(GRAPH_PATTERN* pat, const EDGE_MAP& emap, const int& minsup, pat_fam<GRAPH_PATTERN >& freq_pats, 
                    count_support<GRAPH_PROP, V_Fk1_MINE_PROP, PAT_ST, CC, ALLOC, SM_TYPE>& cs) {
#ifdef PRINT
  cout<<"Fwd extension for "<<pat<<endl;
#endif
//  cout<<"Fwd extension for "<<pat<<endl;
  typedef pat_fam<GRAPH_PATTERN> PAT_FAM;

  const typename GRAPH_PATTERN::RMP_T& rmp=pat->rmost_path();
	pat->random_rmost_path(rmp);
//  const typename GRAPH_PATTERN::RMP_T& rmp=pat->random_rmost_path(); //Random rightmost path
  typename GRAPH_PATTERN::RMP_T::const_iterator rmp_it;
  
  typename EDGE_MAP::CONST_NIT nit;
  typename EDGE_MAP::CONST_LIT lit;  
  GRAPH_PATTERN* edge=0;
  GRAPH_PATTERN* cand_pat=0;
  int lvid;
//Map for storing the can code strings
map<string,int> can_strings;
map<string,int>::iterator can_string_it;
 
  for(rmp_it=rmp.end()-1; rmp_it>=rmp.begin(); rmp_it--) {
//cout<<"rmp_for loop"<<endl;

//The position at which the extension occurs is random..select that randomly
    const typename GRAPH_PATTERN::VERTEX_T& src_v=pat->label(*rmp_it);

//	cout<<"The src label is "<<src_v<<endl;
	if(!emap.has_neighbor(src_v))
	{
//		cout<<"has no neighbor "<<endl;
		continue;
	}
    const typename EDGE_MAP::NEIGHBORS& nbrs=emap.get_neighbors(src_v);
//cout<<"Number of neighbors is "<<nbrs.size();
//randomize the neighbors...
//cout<<"Neighbors before randomization "<<endl;
vector<const char*, ALLOC<const char*> > nbrs_keys;
    for(nit=nbrs.begin(); nit!=nbrs.end(); nit++) {
	nbrs_keys.push_back(nit->first);
//	cout<<"Neighbor is "<<nit->first<<endl;
	}
//Now randomize the keys
random_shuffle(nbrs_keys.begin(),nbrs_keys.end());
//cout<<"Neighbors after randomization "<<endl;
//for(int rand_nbr=0; rand_nbr<nbrs_keys.size(); rand_nbr++) {
//cout<<"Neighbor is "<<nbrs_keys[rand_nbr]<<endl;
//}
    //for(nit=nbrs.begin(); nit!=nbrs.end(); nit++) {
    for(int rand_nbr=0; rand_nbr<nbrs_keys.size(); rand_nbr++) {
//cout<<"nbrs_for loop"<<endl;
	nit = nbrs.find(nbrs_keys[rand_nbr]);
      const typename GRAPH_PATTERN::VERTEX_T& dest_v=nit->first;
//	cout<<"The dest label is "<<dest_v<<endl;

      for(lit=nit->second.begin(); lit!=nit->second.end(); lit++) {
//cout<<"lbls_for loop"<<endl;
        // construct new candidate
//	cout<<"Before adding the edge"<<endl;
        cand_pat=pat->clone();
        lvid=cand_pat->add_vertex(dest_v);
        cand_pat->add_out_edge(*rmp_it, lvid, *lit);
        cand_pat->add_out_edge(lvid, *rmp_it, *lit);
//	cout<<"After adding the edge"<<endl;
#ifdef PRINT
        // cout<<"New candidate="<<cand_pat<<endl;
#endif
//	cout<<"Before isomorphism"<<endl;
        if(!check_isomorphism(cand_pat)) {
#ifdef PRINT
          cout<<"New candidate="<<cand_pat<<endl;
          cout<<"candidate did NOT pass isomorphism"<<endl;
#endif
          delete cand_pat;
          cand_pat=0;
          continue;
        } 
        else {
#ifdef PRINT
          cout<<"New candidate="<<cand_pat<<endl;
          cout<<"candidate passed isomorphism"<<endl;
#endif
        }
//	cout<<"After isomorphism"<<endl;
//Jan 31 2011 : Check if the new pattern constructed is not same as one that was already constructed
//For this : Generate a single huge string from the canonical code and store in a map that is initialized in this function...if the pattern is string is already present in this vector..ignore this pattern...

	string can_string = cand_pat->get_can_string();
	can_string_it = can_strings.find(can_string);
	//cout<<"The candidate string is \n"<<can_string<<endl;
	if(can_string_it != can_strings.end())
	{

	//	cout<<"Pattern string is already present"<<endl;
		//can_string is not present already..
		delete cand_pat;
		cand_pat=0;
		continue;
	}
	else
	{
	//	cout<<"First time encountering this pattern "<<endl;
		can_strings[can_string] = 1;
	}

//
        // create edge
        edge=new GRAPH_PATTERN;
        //if(src_v<dest_v)
//	cout<<"Source edge is "<<src_v<<" and the dest vertex is "<<dest_v<<endl;
          make_edge(edge, src_v, dest_v, *lit);
        //else
          //make_edge(edge, dest_v, src_v, *lit);

#ifdef PRINT
        cout<<"Trying with edge="<<edge<<endl;
#endif
  //      cout<<"Trying with edge="<<edge<<endl;
        // count support
//	cout<<"Before count support"<<endl;
        cs.count(pat, edge, &cand_pat, minsup, 1);
//	cout<<"After count support"<<endl;
     
        delete edge;
        edge=0;
      
        // recursive call for frequent graph
        if(cand_pat->is_valid(minsup)) {
          // freq_pats.push_back(cand_pat);
		//If the candidate is frequent then continue extending..so return from here a signal that continues to the next level..
//          if(print)
 //           cout << cand_pat;
          freq_pats_count++;
          cand_gen(cand_pat, emap, minsup, freq_pats, cs);
          cs.delete_vat(cand_pat);
	//a valid extension..so exit from here with a value of 1
		return 1;
        }
        // else {
//	cout<<"the candidate is not frequent "<<endl;
        delete cand_pat;
        cand_pat=0;
        // }

      }//end for lit
    }//end for nit      
  }//end for rmp_it
//Reaches here if it cannot be extended further...send the message to the caller function...if the sister function back_extension also returns the same signal it means that the pattern cannot be extended on the right most path which means that execution reached one maximal pattern
//since it cannot be extended...return 0
return 0;
}//end fwd_extensions()


/** Populates p with a single-edged pattern;
    v1 is first vertex, v2 is second; 
    It also populates p's canonical code */
// it is assumed that v1<=v2 for VAT intersection to work, this function  does 
// not verify it;
// in particular, storage_manager for a single undirected edge (A-B) stores it 
// as canonical code A-B and not B-A.

template<typename pattern, typename V_T, typename E_T>
void make_edge(pattern* p, const V_T& v1,const V_T& v2, const E_T& e) {
  p->add_vertex(v1);
  p->add_vertex(v2);
  p->add_out_edge(0, 1, e);
  p->add_out_edge(1, 0, e);
  p->init_canonical_code(five_tuple<V_T, E_T>(0, 1, p->label(0), e, p->label(1)));
} //end make_edge()

#endif
