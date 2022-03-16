#ifndef CHITECH_BYTE_ARRAY_H
#define CHITECH_BYTE_ARRAY_H

#include <cstddef>
#include <vector>
#include <stdexcept>

namespace chi_data_types
{
  class ByteArray
  {
  protected:
    std::vector<std::byte> m_raw_data;
    size_t                 m_offset=0;

  public:
    ByteArray() = default;
    explicit
    ByteArray(const size_t raw_data_size) : m_raw_data(raw_data_size) {}
    explicit
    ByteArray(std::vector<std::byte>&& raw_data) : m_raw_data(std::move(raw_data)) {}
    explicit
    ByteArray(const std::vector<std::byte>& raw_data) : m_raw_data(raw_data) {}


    /**Uses the template type T to convert an associated value (of type T)
     * to a sub-array of std::bytes and adds it to the internal byte-array.
     *
     * The template type T must support sizeof.*/
    template<typename T> void Write(const T& value)
    {
      const size_t num_bytes = sizeof(T);
      const std::byte* value_byte_array =
        reinterpret_cast<const std::byte*>(&value);

      for (size_t b=0; b<num_bytes; ++b)
        m_raw_data.push_back(value_byte_array[b]);
    }

    /**Uses the template type `T` to convert `sizeof(T)` number of bytes to
     * a value of type `T`. The bytes are pulled from the internal byte-array
     * starting at the internal address `m_offset` which can be viewed with a
     * call to `Offset()`. This offset is incremented by the amount of bytes
     * used and can be reset with a call to `Seek(0)`.
     *
     * Bounds-checking is performed by checking if the internal byte array
     * has the required number of bytes available. If this check fails then
     * this call will return an `out_of_range` exception.*/
    template<typename T> T Read()
    {
      const size_t num_bytes = sizeof(T);
      if ((m_offset + num_bytes - 1) >= m_raw_data.size())
        throw std::out_of_range(
          std::string("ByteArray reading error. ") +
        " Typename: " + std::string(typeid(T).name()) +
        " m_offset: " + std::to_string(m_offset) +
        " size: " + std::to_string(m_raw_data.size()) +
        " num_bytes to read: " + std::to_string(num_bytes));

      T value = *reinterpret_cast<T*>(&m_raw_data[m_offset]);
      m_offset += num_bytes;

      return value;
    }

    /**Uses the template type T to convert sizeof(T) number of bytes to
     * a value of type T. The bytes are pulled from the internal byte-array
     * starting at the internal address specified by the argument "address".
     * An optional argument next_address can be used to return the location
     * of the value after the value read.
     * This offset is the given "address" argument incremented by the amount
     * of bytes used.
     *
     * Bounds-checking is performed by checking if the internal byte array
     * has the required number of bytes available. If this check fails then
     * this call will return an `out_of_range` exception.*/
    template<typename T> T Read(const size_t address,
                                size_t* next_address = nullptr) const
    {
      const size_t num_bytes = sizeof(T);
      if ((address + num_bytes - 1) >= m_raw_data.size())
        throw std::logic_error(
          std::string("ByteArray reading error. ") +
        " Typename: " + std::string(typeid(T).name()) +
        " address: " + std::to_string(address) +
        " size: " + std::to_string(m_raw_data.size()) +
        " num_bytes to read: " + std::to_string(num_bytes));

      T value = *reinterpret_cast<const T*>(&m_raw_data[address]);
      if (next_address != nullptr) *next_address = address + num_bytes;

      return value;
    }

    /**Appends a `ByteArray` to the current internal byte array.*/
    void Append(const ByteArray& other_raw)
    {
      auto& master = m_raw_data;

      const auto& slave = other_raw.Data();
      for (std::byte slave_value : slave)
        master.push_back(slave_value);
    }

    /**Appends bytes from a `std::vector<std::byte>` to the internal
     * byte array.*/
    void Append(const std::vector<std::byte>& other_raw)
    {
      auto& master = m_raw_data;

      const auto& slave = other_raw;
      for (std::byte slave_value : slave)
        master.push_back(slave_value);
    }

    /**Clears the internal byte array and resets the address offset to zero.*/
    void Clear()
    {
      m_raw_data.clear();
      m_offset = 0;
    }

    /**Moves the address marker to the supplied address.*/
    void Seek(const size_t address=0) {m_offset = address;}
    /**Returns the internal address marker position.*/
    size_t Offset() const {return m_offset;}
    /**Determines if the internal address marker is beyond the internal
     * byte array.*/
    bool EndOfBuffer() const {return m_offset >= m_raw_data.size();}
    /**Returns the current size of the internal byte array.*/
    size_t Size() const {return m_raw_data.size();}

    /**Returns a reference to the internal byte array.*/
    std::vector<std::byte>& Data() {return m_raw_data;}
    /**Returns a const reference of the internal byte array.*/
    const std::vector<std::byte>& Data() const {return m_raw_data;}
  };
}//namespace chi_data_types


#endif //CHITECH_BYTE_ARRAY_H
