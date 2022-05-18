#ifndef CHITECH_CHI_LOG_EXCEPTIONS_H
#define CHITECH_CHI_LOG_EXCEPTIONS_H

#define ChiInvalidArgument(condition, message)  \
{                                               \
  if (condition)                                \
    throw std::invalid_argument(                \
      std::string(__PRETTY_FUNCTION__) + ": " + \
      #message);                                \
}

#define ChiLogicalError(condition, message)     \
{                                               \
  if (condition)                                \
    throw std::invalid_argument(                \
      std::string(__PRETTY_FUNCTION__) + ": " + \
      #message);                                \
}

#endif //CHITECH_CHI_LOG_EXCEPTIONS_H
