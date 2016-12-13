// by Wei Shi (based on code by Mark Weiss)
// last modified 09/23/16
// Modified DisjointSets data structure for use in image processing

#ifndef DISJ_SETS_H
#define DISJ_SETS_H

// DisjSets class
//
// CONSTRUCTION: with int representing initial number of sets
//
// ******************PUBLIC OPERATIONS*********************
// void union( root1, root2 ) --> Merge two sets
// int find( x )              --> Return set containing x
// ******************ERRORS********************************
// No error checking is performed
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>
#include <unordered_set>

using namespace std;

/**
 * Disjoint set class.
 * Use union by rank and path compression.
 * Elements in the set are numbered starting at 0.
 */

class DisjointSets
{
  public:
		DisjointSets();
    explicit DisjointSets( size_t numElements );
    //explicit DisjointSets( unsigned short numElements );

		// looks for the root of x
    int find( const int& x ) const;
 		int find( int x );

		size_t size() const;

		// @param index is the index of the element to check if it is in use
		bool isUsed(size_t index) const;

		// always use unionSets(find(x), find(y)) to find the root of x and root of y
		// before unioning two sets.
    void unionSets( int root1, int root2 );

		void print() const;
		void printGroup() const;
		void printValidIndex() const;

		// post-condition sets valid_index vector
		void setValidIndex(vector<bool>& rhs);

		const vector<int>& getDisjointSet() const;
		DisjointSets& operator=(const DisjointSets& rhs);
		
  private:
		vector<bool> valid_index;
    vector<int> s;
};

#endif
