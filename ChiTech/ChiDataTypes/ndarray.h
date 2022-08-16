#ifndef CHITECH_NDARRAY_H
#define CHITECH_NDARRAY_H

#include <cstddef>
#include <type_traits>
#include <cassert>
#include <vector>
#include <array>
#include <stdexcept>

namespace chi_data_types
{

template<typename T>
class NDArray
{
private:
  size_t                m_rank;
  size_t*               m_dimensions;
  size_t*               m_strides;
  size_t                m_size;
  T*                    m_base;

  template<bool...>
  struct bool_pack{};

  template<class... U>
  using conjunction = std::is_same<
                      bool_pack<true,U::value...>,
                      bool_pack<U::value..., true>
                      >;

  template<typename... U>
  using AllIntegral = typename conjunction<std::is_integral<U>...>::type;

public:
  /** Creates an array with the specified number of elements in each dimension,
   *  from a vector-list.
   *
   *  \param args `std::vector` list of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template<typename D>
  explicit
  NDArray(const std::vector<D>& dims) :
    m_rank(0), m_dimensions(nullptr), m_strides(nullptr),
    m_size(0), m_base(nullptr)
  {
    static_assert(std::is_integral<D>::value,
      "NDArray dims argument must have vector of integral types.");
    const size_t N = dims.size();

    m_rank = N;

    m_strides = new size_t[N];
    m_dimensions = new size_t[N];

    //Populate dimensions
    for(size_t i = 0; i < N; ++i)
      m_dimensions[i] = dims[i];

    //Populate strides and size
    m_size = 1;
    for(size_t i = 0; i < N; ++i)
    {
      assert(m_dimensions[i] > 0);
      m_size *= m_dimensions[i];
      m_strides[i] = 1;
      for(size_t j = i+1; j < N; ++j)
        m_strides[i] *= m_dimensions[j];
    }

    m_base = new T[m_size];
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from an array.
   *
   *  \param args `std::array` list of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template<typename D, size_t N>
  explicit
  NDArray(const std::array<D, N>& dims) :
    m_rank(0), m_dimensions(nullptr), m_strides(nullptr),
    m_size(0), m_base(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have array of integral types.");
    m_rank = N;

    m_strides = new size_t[N];
    m_dimensions = new size_t[N];

    //Populate dimensions
    for(size_t i = 0; i < N; ++i)
      m_dimensions[i] = dims[i];

    //Populate strides and size
    m_size = 1;
    for(size_t i = 0; i < N; ++i)
    {
      assert(m_dimensions[i] > 0);
      m_size *= m_dimensions[i];
      m_strides[i] = 1;
      for(size_t j = i+1; j < N; ++j)
        m_strides[i] *= m_dimensions[j];
    }

    m_base = new T[m_size];
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from an initializer-list.
   *
   *  \param args `std::vector` list of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template<typename D>
  NDArray(const std::initializer_list<D>& dims) :
    m_rank(0), m_dimensions(nullptr), m_strides(nullptr),
    m_size(0), m_base(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");
    const size_t N = dims.size();

    m_rank = N;

    m_strides = new size_t[N];
    m_dimensions = new size_t[N];

    //Populate dimensions
    {
      size_t i=0;
      for (size_t val : dims)
      {
        m_dimensions[i] = val; ++i;
        if (i >= N) break;
      }
    }

    //Populate strides and size
    m_size = 1;
    for(size_t i = 0; i < N; ++i)
    {
      assert(m_dimensions[i] > 0);
      m_size *= m_dimensions[i];
      m_strides[i] = 1;
      for(size_t j = i+1; j < N; ++j)
        m_strides[i] *= m_dimensions[j];
    }

    m_base = new T[m_size];
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from a vector. Each entry in the array is assigned the designated value.
   *
   *  \param args `std::vector` list of the number of elements in each
   *              dimension.
   *  \param value The value to assing to each element.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template<typename D>
  explicit
  NDArray(const std::vector<D>& dims, T value) :
    m_rank(0), m_dimensions(nullptr), m_strides(nullptr),
    m_size(0), m_base(nullptr)
  {
    static_assert(std::is_integral<D>::value,
      "NDArray dims argument must have vector of integral types.");
    const size_t N = dims.size();

    m_rank = N;

    m_strides = new size_t[N];
    m_dimensions = new size_t[N];

    //Populate dimensions
    for(size_t i = 0; i < N; ++i)
      m_dimensions[i] = dims[i];

    //Populate strides and size
    m_size = 1;
    for(size_t i = 0; i < N; ++i)
    {
      assert(m_dimensions[i] > 0);
      m_size *= m_dimensions[i];
      m_strides[i] = 1;
      for(size_t j = i+1; j < N; ++j)
        m_strides[i] *= m_dimensions[j];
    }

    m_base = new T[m_size];
    for (size_t i=0; i<m_size; ++i)
      m_base[i] = value;
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from an array. Each entry in the array is assigned the designated value.
   *
   *  \param args `std::array` list of the number of elements in each
   *              dimension.
   *  \param value The value to assing to each element.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template<typename D, size_t N>
  explicit
  NDArray(const std::array<D, N>& dims, T value) :
    m_rank(0), m_dimensions(nullptr), m_strides(nullptr),
    m_size(0), m_base(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");
    m_rank = N;

    m_strides = new size_t[N];
    m_dimensions = new size_t[N];

    //Populate dimensions
    for(size_t i = 0; i < N; ++i)
      m_dimensions[i] = dims[i];

    //Populate strides and size
    m_size = 1;
    for(size_t i = 0; i < N; ++i)
    {
      assert(m_dimensions[i] > 0);
      m_size *= m_dimensions[i];
      m_strides[i] = 1;
      for(size_t j = i+1; j < N; ++j)
        m_strides[i] *= m_dimensions[j];
    }

    m_base = new T[m_size];
    for (size_t i=0; i<m_size; ++i)
      m_base[i] = value;
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from an initializer-list. Each entry in the array is assigned the
   *  designated value.
   *
   *  \param args `std::initializer` list of the number of elements in each
   *              dimension.
   *  \param value The value to assing to each element.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template<typename D>
  NDArray(const std::initializer_list<D>& dims, T value) :
    m_rank(0), m_dimensions(nullptr), m_strides(nullptr),
    m_size(0), m_base(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");
    const size_t N = dims.size();

    m_rank = N;

    m_strides = new size_t[N];
    m_dimensions = new size_t[N];

    //Populate dimensions
    {
      size_t i=0;
      for (size_t val : dims)
      {
        m_dimensions[i] = val; ++i;
        if (i >= N) break;
      }
    }

    //Populate strides and size
    m_size = 1;
    for(size_t i = 0; i < N; ++i)
    {
      assert(m_dimensions[i] > 0);
      m_size *= m_dimensions[i];
      m_strides[i] = 1;
      for(size_t j = i+1; j < N; ++j)
        m_strides[i] *= m_dimensions[j];
    }

    m_base = new T[m_size];
    for (size_t i=0; i<m_size; ++i)
      m_base[i] = value;
  }

  /** Creates an empty array.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an empty array and initializes the reference count to one.
   */
  NDArray() :
    m_rank(0), m_dimensions(nullptr), m_strides(nullptr),
    m_size(0), m_base(nullptr)
  {
  }

  /** Copy construct from another array.
   *  \param other The array to copy.
   */
  NDArray(NDArray<T> const &other) :
    m_rank(other.m_rank),
    m_dimensions(nullptr),
    m_strides(nullptr),
    m_size(other.m_size),
    m_base(nullptr)
  {

    const size_t N = m_rank;
    m_strides = new size_t[N];
    m_dimensions = new size_t[N];
    m_base = new T[m_size];

    for(size_t i = 0; i < N; ++i)
    {
      m_dimensions[i] = other.m_dimensions[i];
      m_strides[i]    = other.m_strides[i];
    }

    for (size_t i=0; i<m_size; ++i)
      m_base[i] = other.m_base[i];
  }

  /** Assign from another array.
   *  \param other The array to copy.
   */
  NDArray<T> & operator=(NDArray<T> const &other)
  {
    NDArray<T>(other).swap(*this);
    return *this;
  }

  /**Move constructor.*/
  NDArray(NDArray<T>&& other) noexcept :
    m_rank(std::move(other.m_rank)),
    m_dimensions(std::move(other.m_dimensions)),
    m_strides(std::move(other.m_strides)),
    m_size(std::move(other.m_size)),
    m_base(std::move(other.m_base))
  {
  }

  /**Deleted move assignment operator*/
  NDArray<T> & operator=(NDArray<T> &&) = delete;

  /**Sets a value to all the items in the array.*/
  void set(T value)
  {
    for (size_t i = 0; i < m_size; ++i)
      m_base[i] = value;
  }

  /** Swap the contents of this array with another array.
   *  \param other The array to swap with.
   */
  void swap(NDArray<T> &other)
  {
    std::swap(m_rank, other.m_rank);
    std::swap(m_dimensions, other.m_dimensions);
    std::swap(m_strides, other.m_strides);
    std::swap(m_size, other.m_size);
    std::swap(m_base, other.m_base);
  }

  /** Returns the number of elements in the array.*/
  size_t size() const noexcept { return m_size; }

  /** Returns true if the array has no elements.*/
  bool empty() const noexcept { return m_size == 0; }

  /** Returns an iterator pointing to the beginning of the array.*/
  T * begin() const noexcept { return m_base; }

  /** Returns a constant iterator pointing to the beginning of the array.*/
  const T * cbegin() const noexcept { return m_base; }

  /** Returns an iterator pointing to the end of the array.*/
  T * end() const noexcept { return m_base + m_size; }

  /** Returns a constant iterator pointing to the end of the array*/
  const T * cend() const noexcept { return m_base + m_size; }

  /** Returns a pointer to the underlying array data.*/
  T * data() const noexcept { return m_base; }

  /** Returns the rank of the array.*/
  size_t rank() const noexcept { return m_rank; }

  /**Returns the dimension of the array.*/
  std::vector<size_t> dimension() const
  {
    std::vector<size_t> dim(m_rank, 0);
    for (size_t i=0; i<m_rank; ++i)
      dim[i] = m_dimensions[i];

    return dim;
  }

  /** Resizes the array with a vector.
   *
   *  \param args std::vector of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This method resizes the array to the specified number of elements. If the
   *  current size is equal to the new size, no memory allocation occurs.
   */
  template<typename D>
  void resize(const std::vector<D>& dims)
  {
    static_assert(std::is_integral<D>::value,
      "NDArray dims argument must have vector of integral types.");

    delete [] m_dimensions;
    m_dimensions = nullptr;
    delete [] m_strides;
    m_strides = nullptr;
    delete [] m_base;
    m_base = nullptr;

    NDArray<T>(dims).swap(*this);
  }

  /** Resizes the array with an array.
   *
   *  \param args std::array of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This method resizes the array to the specified number of elements. If the
   *  current size is equal to the new size, no memory allocation occurs.
   */
  template<typename D, size_t N>
  void resize(const std::array<D, N>& dims)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");

    delete [] m_dimensions;
    m_dimensions = nullptr;
    delete [] m_strides;
    m_strides = nullptr;
    delete [] m_base;
    m_base = nullptr;

    NDArray<T>(dims).swap(*this);
  }

  /** Resizes the array with an initializer_list.
   *
   *  \param args std::initializer_list of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This method resizes the array to the specified number of elements. If the
   *  current size is equal to the new size, no memory allocation occurs.
   */
  template<typename D>
  void resize(const std::initializer_list<D>& dims)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");

    delete [] m_dimensions;
    m_dimensions = nullptr;
    delete [] m_strides;
    m_strides = nullptr;
    delete [] m_base;
    m_base = nullptr;

    NDArray<T>(dims).swap(*this);
  }

  /** Accesses the specified element.
   *  \param args The indices of the desired element.
   *  \return Read/write reference to the element.
   */
  template<typename... Args>
  T& operator()(Args... args) noexcept
  {
    static_assert(AllIntegral<Args...>::value,
      "NDArray::operator[]: All parameters must be of integral type");

    size_t indices[] { static_cast<size_t>(args)... };

    const size_t N = m_rank;

    T *address = m_base + indices[N-1];
    for(size_t i = 0;i < N-1;++i)
      address += m_strides[i] * indices[i];
    return *(address);
  }

  /** Accesses the specified element.
   *  \param args The indices of the desired element.
   *  \return Read reference to the element.
   */
  template<typename... Args>
  T const& operator()(Args... args) const noexcept
  {
    static_assert(AllIntegral<Args...>::value,
      "NDArray::operator[]: All parameters must be of integral type");

    size_t indices[] { static_cast<size_t>(args)... };

    const size_t N = m_rank;

    T *address = m_base + indices[N-1];
    for(size_t i = 0;i < N-1;++i)
      address += m_strides[i] * indices[i];
    return *(address);
  }

  /** Accesses the specified element with safety checks.
   *
   *  \param args The indices of the desired element.
   *  \throw std::invalid_argument if the number of arguments are incorrect and
   *  std::out_of_range if one of the dimension-indices are out of range.
   *  \return Read/write reference to the element.
   */
  template<typename... Args>
  T& at(Args... args) noexcept
  {
    static_assert(AllIntegral<Args...>::value,
      "NDArray::at(): All parameters must be of integral type");

    if (sizeof...(args) != m_rank)
      throw std::invalid_argument("NDArray::at(): Number of arguments " +
        std::to_string(sizeof...(args)) + " not equal to rank " +
        std::to_string(m_rank));

    const size_t N = m_rank;
    size_t indices[] { static_cast<size_t>(args)... };
    for (size_t i=0; i<N; ++i)
      if (indices[i] >= m_dimensions[i])
        throw std::out_of_range("NDArray::at(): Index " + std::to_string(i) +
          " out of range " + std::to_string(indices[i]) +
          " must be <" + std::to_string(m_dimensions[i]));

    T *address = m_base + indices[N-1];
    for(size_t i = 0;i < N-1;++i)
      address += m_strides[i] * indices[i];
    return *(address);
  }

  /** Deletes the array.
   *
   *  The destructor deletes the underlying array data.
   */
  ~NDArray()
  {
    delete [] m_dimensions;
    m_dimensions = nullptr;
    delete [] m_strides;
    m_strides = nullptr;
    delete [] m_base;
    m_base = nullptr;
  }
};

}//namespace chi_data_types

#endif //CHITECH_NDARRAY_H
