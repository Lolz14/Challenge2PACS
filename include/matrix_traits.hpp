#ifndef MATRIX_TYPE_TRAITS
#define MATRIX_TYPE_TRAITS

#include <vector>
#include <array>
#include <map>
#include <type_traits>
#include <concepts>

/*!
*Types and structures definitions to Sparse Matrix.
*/


/*!
* Enumerator that specifies ordering.
*/
enum StorageOrder{
    RowWise = 0,
    ColumnWise = 1
};


/*!
* Enumerator that specifies norms.
*/
enum NormType{
    One = 0,
    Infinity = 1,
    Frobenius = 2
};


/*!
* Tag dispatching implementation for a more efficient handling of NormType selection.
* @param N: template parameters mapping the correct norm as defined in NormType
*/
template <NormType N>
using NT = std::integral_constant<NormType,N>;


/*!
* Definition of the key containing the map for the uncompressed storage.
*/
using key_type = std::array<std::size_t,2>;


/*!
* Functor that defines ColumnWise ordering: to be passed as template parameter.
* @param lhs: const l-value reference to an element of type key_type
* @param rhs: const l-value reference to an element of type key_type
* @return true by comparing column number first, then by comparing row number if column nunmbers are equal
*/
struct IsLessColWise{
    bool operator() (const key_type &lhs, const key_type &rhs) const{
        return ((lhs[1]<rhs[1]) or ((lhs[1]==rhs[1]) and (lhs[0]<rhs[0]))); 
    };
};


/*!
* Template definition for  uncompressed matrix.
* @param T: the type of elements that are stored
* @param O: StorageOrder 
* @return using std::conditional: O  is 0 (RowWise): the map will have the defaulted lexicographical less operator;
*                                 O is 1 (ColumnWise): the map will have the IsLessColWise functor as less operator
*/
template <class T,StorageOrder O>
using MatrixUncompressed = std::conditional<O,
                                            std::map<key_type, T, IsLessColWise>,
                                            std::map<key_type,T>>::type;


#endif ///MATRIX_TYPE_TRAITS