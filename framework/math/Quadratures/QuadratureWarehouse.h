#ifndef CHITECH_QUADRATUREWAREHOUSE_H
#define CHITECH_QUADRATUREWAREHOUSE_H

#include "quadrature.h"

#include <map>

namespace chi_math
{

/**Singleton object for the easy reuse of quadratures. This singleton
* has been created to remove the need for spatial discretizations to store
* multiple quadratures at construction.*/
class QuadratureWarehouse
{
public:
  typedef std::pair<std::string, QuadratureOrder> QuadratureKey;

  /**Returns the singleton instance of the QuadratureWarehouse*/
  static QuadratureWarehouse& GetInstance() noexcept;

  // Deleted copy constructor
  QuadratureWarehouse(const QuadratureWarehouse&) = delete;
  // Deleted assignment operator
  QuadratureWarehouse operator=(const QuadratureWarehouse&) = delete;

  /**Either builds a new quadrature of the given type or retrieves an
  * existing one.*/
  template <class T>
  const Quadrature&
  GetQuadrature(QuadratureOrder q_order)
  {
    const std::string type_name = typeid(T).name();
    const auto key = std::make_pair(type_name, q_order);
    auto it = quadratures_.find(key);

    if (it != quadratures_.end()) return *it->second;
    else
    {
      auto result =
        quadratures_.insert(std::make_pair(key, std::make_unique<T>(q_order)));

      auto& quad = *result.first->second;
      return quad;
    }
  }

  /**Either builds a new quadrature of the given type or retrieves an
  * existing one. This overload allows for changing the range of 1-dimensional
  * quadrature.*/
  template <class T>
  const Quadrature&
  GetQuadrature(QuadratureOrder q_order,
                const std::pair<double, double> range)
  {
    const std::string type_name = typeid(T).name();
    const auto key = std::make_pair(type_name, q_order);
    auto it = quadratures_.find(key);

    if (it != quadratures_.end()) return *it->second;
    else
    {
      auto result =
        quadratures_.insert(std::make_pair(key, std::make_unique<T>(q_order)));

      auto& quad = *result.first->second;
      quad.SetRange(range);
      return quad;
    }
  }

private:
  QuadratureWarehouse() = default;
  std::map<QuadratureKey, std::unique_ptr<Quadrature>> quadratures_;
};

} // namespace chi_math

#endif // CHITECH_QUADRATUREWAREHOUSE_H
