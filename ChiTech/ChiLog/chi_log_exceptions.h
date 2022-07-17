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
    throw std::logic_error(                     \
      std::string(__PRETTY_FUNCTION__) + ": " + \
      #message);                                \
}

#define ChiRecoverableInvalidArgument(condition, message)  \
{                                                    \
  if (condition)                                     \
    throw std::RecoverableException(                 \
      std::string("Recoverable Invalid Argument: "), \
      std::string(__PRETTY_FUNCTION__) + ": " +      \
      #message);                                     \
}

#define ChiRecoverableLogicalError(condition, message)     \
{                                               \
  if (condition)                                \
    throw std::RecoverableException(            \
      std::string("Recoverable Logic Error: ")  \
      std::string(__PRETTY_FUNCTION__) + ": " + \
      #message);                                \
}

#endif //CHITECH_CHI_LOG_EXCEPTIONS_H
