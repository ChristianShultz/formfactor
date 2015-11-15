#ifndef ERR_MAC_H_H_GUARD
#define ERR_MAC_H_H_GUARD

namespace splashMessage
{

  namespace splash
  {
    void splash(const char *prettyFunction, const char *file, const int line, const char *msg);
  }
}

#define SPLASH(msg) \
  do \
    { \
      splashMessage::splash::splash(__PRETTY_FUNCTION__, __FILE__, __LINE__, (msg) ); \
    } while (0)

#endif
