#ifndef TOKENIZE_H
#define TOKENIZE_H 

#include <vector> 
#include <string>
#include <sstream>


namespace radmat
{

  template<typename T> 
    std::vector<std::string> tokenize(const std::string &s,
        T &delim_t ,
        const bool keep_empty_tokens=false)
    {
      std::vector<std::string> tokens;

      std::stringstream ss; 
      ss << delim_t; 

      std::string delim = ss.str();

      if(delim.empty())
        return std::vector<std::string>(1,s);

      std::string::const_iterator front,back;
      front = s.begin();
      const std::string::const_iterator end = s.end();

      while(true)
      {

        back = std::search(front,end,delim.begin(),delim.end());
        std::string token(front,back);
        if(keep_empty_tokens || !!!token.empty())
          tokens.push_back(token);

        if(back == end)
          break;

        front = back + delim.size();
      }

      return tokens; 
    }

} // radmat

#endif /* TOKENIZE_H */
