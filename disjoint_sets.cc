// by Wei Shi (based on code by Mark Weiss)
// last modified 09/23/16
// Modified DisjointSets data structure for use in image processing

#include "disjoint_sets.h"


DisjointSets::DisjointSets() {
   
}

/**
 * Construct the disjoint sets object.
 * numElements is the initial number of disjoint sets.
 */

DisjointSets::DisjointSets( size_t numElements ) : s( numElements ) {
    for( int i = 0; i < s.size( ); i++ )
        s[ i ] = -1;
}

/**
 * Union two disjoint sets.
 * For simplicity, we assume root1 and root2 are distinct
 * and represent set names.
 * root1 is the root of set 1.
 * root2 is the root of set 2.
 */

void DisjointSets::unionSets( int root1, int root2 ) {

	// check to see if they are already in the same set:
	if(root1 == root2) {
		return;

	}
	// e.g., root1 = 3, root2 = 6
  if( s[ root2 ] < s[ root1 ] ) {  // root2 is deeper
		// Make root2 new root
    s[ root1 ] = root2;   


	}

  else
  {
    if( s[ root1 ] == s[ root2 ] )
        s[ root1 ]--;          // Update height if same

		// s[6] = 3
    s[ root2 ] = root1;        // Make root1 new root of root2
  }
}

/**
 * Perform a find.
 * Error checks omitted again for simplicity.
 * Return the set containing x.
 */

int DisjointSets::find( const int& x ) const {
    if( s[ x ] < 0 )
        return x;
    else
        return find( s[ x ] );
}


/**
 * Perform a find with path compression.
 * Error checks omitted again for simplicity.
 * Return the set containing x.
 */

int DisjointSets::find( int x ) {


		// x = 3 s[x]
    if( s[ x ] < 0 )
        return x;
    else
        return s[ x ] = find( s[ x ] );
}

size_t DisjointSets::size() const {

	return s.size();


}

bool DisjointSets::isUsed(size_t index) const{

	return valid_index[index];

}

void DisjointSets::print() const {
	for(auto& i: s) {
		cout << i << ' ';
	}
	cout << endl;

}

void DisjointSets::setValidIndex(vector<bool>& rhs) {

	valid_index = rhs;
	
}

void DisjointSets::printValidIndex() const{

	for(const bool& i: valid_index) {

		cout << i << ' ';
	}

	cout << endl;
}

const vector<int>& DisjointSets::getDisjointSet() const{
	return s;
}
DisjointSets& DisjointSets::operator=(const DisjointSets& rhs) {

	valid_index = rhs.valid_index;
	s = rhs.s;
	//group = rhs.group;



}
















