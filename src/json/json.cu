#include "json.h"

void json::print(std::ostream& out, Value& obj)
{
  struct V : json::Value::Visitor
  {
    V(std::ostream& o) : out(o), tab_level(0) { }

    void visit(const json::Object& obj)
    {
      out << "{" << std::endl;
      ++tab_level;
      for(json::Object::const_iterator it = obj.begin(); it != obj.end(); ++it)
      {
        for(int i = 0; i < tab_level; ++i)
          out << "\t";
        out << "\"" << it->first << "\":";
        it->second->accept(*this);
        if(++it != obj.end())
          out << "," << std::endl;
        --it;
      }
      out << std::endl << "}" << std::endl;
      --tab_level;
    }
    void visit(const json::Array& obj)
    {
      out << "[";
      for(json::Array::const_iterator it = obj.begin(); it != obj.end(); ++it)
      {
        (*it)->accept(*this);
        if(++it != obj.end())
          out << ",";
        --it;
      }
      out << "]";
    }
    void visit(const json::String& obj) { out << obj; }
    void visit(const json::Number& obj) { out << obj.val; }
    void visit(const json::Bool& obj) { out << (obj.val?"true":"false"); }
    void visit(const json::Null& obj) { out << "null"; }
    int tab_level;
    std::ostream& out;
  };
  V vis(out);
  obj.accept(vis);
}

json::Object json::parse(std::istream& in, int max_depth)
{
  JSON_checker jc = new_JSON_checker(max_depth);
  char cur;
  int line = 1;
  json::Object o;
  while(in >> cur)
  {
    if(JSON_checker_char(jc,cur))
    {
      if(cur == '\n')
        ++line;

      return *parse_object(in,jc,line);
    }else
    {
      std::cout << "Invalid JSON Syntax. at line: " << line
                << " character: " << cur << std::endl;
      exit(1);
    }
  }
  exit(1);
  return o;
}

json::Object* json::parse_object(std::istream& in, JSON_checker& jc, int& line)
{
  char cur;
  json::Object* document = new json::Object();
  std::string temp_key("");
  std::string temp_num("");
  json::Value * temp_value;
  bool have_key = false;
  bool have_value = false;

  while(in >> cur)
  {
    if(JSON_checker_char(jc,cur))
    {
      int state = jc->state;

      switch(state)
      {
        case 1: // OK
          if(temp_num.size() != 0)
          {
            json::Number* num = new json::Number(temp_num);
            temp_value = num;
            temp_num = "";
          }
          have_value = true;
          break;
        case 2: // Object
          temp_value = json::parse_object(in,jc,line);
          break;
        case 3: // Key
          if(temp_num.size() != 0)
          {
            json::Number* num = new json::Number(temp_num);
            temp_value = num;
            temp_num = "";
            have_value = true;
          }
          break;
        case 4: // Colon
          have_key = true;
          break;
        case 6: // Array
          temp_value = json::parse_array(in,jc,line);
          have_value = true;
          break;
        case 7: // String
          if(cur != '\"')
            temp_key += cur;
          break;
        case 13: case 14: case 15:
        case 16: case 17: case 18:
        case 19: // Number
          temp_num += cur;
          break;
        case 22: // True
        {
          json::Bool* b_1 = new json::Bool(true);
          temp_value = b_1;
          break;
        }
        case 26: // False
        {
          json::Bool* b2 = new json::Bool(false);
          temp_value = b2;
          break;
        }
        case 29: // Null
          json::Null* n = new json::Null();
          temp_value = n;
      }

      if(have_key and have_value)
      {
        (*document)[temp_key] = temp_value;
        temp_key = "";
        have_key = false;
        have_value = false;
      }

      if(cur == '}')
        break;
    }else
    {
      std::cout << "Invalid JSON Syntax. at line: " << line
                << " character: " << cur << std::endl;
      exit(1);
    }
  }
  return document;
}

json::Array* json::parse_array(std::istream& in, JSON_checker& jc, int& line)
{
  char cur;
  json::Array* array = new json::Array();
  std::string temp_num("");
  std::string temp_string("");
  json::Value* temp_value;
  bool have_value = false;

  while(in >> cur)
  {
    if(JSON_checker_char(jc,cur))
    {
      int state = jc->state;

      switch(state)
      {
        case 1: // OK
          if(temp_num.size() != 0)
          {
            json::Number* num = new json::Number(temp_num);
            temp_value = num;
            temp_num = "";
          }
          have_value = true;
          break;
        case 2: // Object
          temp_value = json::parse_object(in,jc,line);
          break;
        case 5: // Value
          if(temp_num.size() != 0)
          {
            json::Number* num = new json::Number(temp_num);
            temp_value = num;
            temp_num = "";
            have_value = true;
          }
          break;
        case 6: // Array
          temp_value = json::parse_array(in,jc,line);
          have_value = true;
          break;
        case 7: // String
          if(cur != '\"')
            temp_string += cur;
          break;
        case 13: case 14: case 15:
        case 16: case 17: case 18:
        case 19: // Number
          temp_num += cur;
          break;
        case 22: // True
        {
          json::Bool* b_1 = new json::Bool(true);
          temp_value = b_1;
          break;
        }
        case 26: // False
        {
          json::Bool* b2 = new json::Bool(false);
          temp_value = b2;
          break;
        }
        case 29: // Null
          json::Null* n = new json::Null();
          temp_value = n;
      }

      if(have_value)
      {
        array->push_back(temp_value);
        have_value = false;
      }

      if(cur == ']')
        break;
    }else
    {
      std::cout << "Invalid JSON Syntax. at line: " << line
                << " character: " << cur << std::endl;
      exit(1);
    }
  }
  return array;
}
