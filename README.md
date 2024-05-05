# Sparse Matrix Library

Welcome to the Sparse Matrix Library, a C++ toolkit crafted to manage sparse matrices with efficiency and flexibility. Whether you're handling large datasets or intricate computational problems, this library equips you with the tools needed to streamline your operations.

## Data Structure

The core class in this library is `Matrix<T,O>` under the namespace `algebra`, which is a template class. The template parameter `T` refers to the type of values stored in the matrix (integral, floating-point, or complex types), while the parameter `O` indicates the storage order (`RowWise` or `ColumnWise`). Once a matrix is constructed with a storage order, it cannot be changed.

### Formats

- **Uncompressed Format (COOmap)**: Data stored using `std::map`, with keys as `std::array<std::size_t,2>` representing row and column indices, and values as type `T`. Row-wise storage is handled using lexicographical ordering of the keys, while column-wise storage is handled using a functor.
  
- **Compressed Format (CSR/CSC)**: Data stored using three `std::vector`s: one for values of type `T` and the other two for index storage of type `std::size_t`.

## Requirements


- PACS utility `Chrono.hpp`.
- Modification of `Makefile` required: set the proper `PACS_ROOT` at the beginning.
- Modification of `Doxyfile` required: set `INCLUDE_PATH = ./include \ .\ PACS_ROOT/include`.

## Usage

### Compilation

Compile the code with optimization flags already set to -O3

```
make
```

### Execution

Run the compiled executable:

```
./main
```

### Documentation

Generate documentation using Doxygen:

```
make doc
```

This creates a `doc` directory containing HTML and LaTeX subfolders. Open `index.html` in the `html` folder to view the documentation in a web browser.

## Files

- **main.cpp**: Demonstrates various operations using matrices read from an MMF file and evaluates their performance.
- **/include**: Contains header files:
  - `matrix.hpp`: Declaration of the `Matrix` and definition of some methods.
  - `utils.hpp`: Definition of all the methods not defined in `matrix.hpp`.
  - `csv.h`: Declaration and definition of the csv class.
  - `matrix_traits.hpp`: Type Traits used in the implementation.




