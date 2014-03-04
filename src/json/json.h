#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "JSON_checker.h"

#ifndef JSON_H
#define JSON_H

namespace json
{

struct Value
{
  struct Visitor;

  virtual void accept(Visitor&) const = 0;
  virtual ~Value() { };
};

struct Object;
struct Array;
struct String;
struct Number;
struct Bool;
struct Null;

struct Value::Visitor
{
  //virtual void visit_value(const Value&);
  virtual void visit(const Object&) = 0;
  virtual void visit(const Array&) = 0;
  virtual void visit(const String&) = 0;
  virtual void visit(const Number&) = 0;
  virtual void visit(const Bool&) = 0;
  virtual void visit(const Null&) = 0;
};

struct Object : std::map<std::string,Value*>, Value
{
  virtual void accept(Value::Visitor& v) const
  {
    v.visit(static_cast<const Object&>(*this));
  }

  virtual ~Object()
  {
    for(std::map<std::string,Value*>::iterator it = this->begin();
        it != this->end(); ++it)
    {
      delete it->second;
    }
  }
};

struct Array : std::vector<Value*>, Value
{
  virtual void accept(Value::Visitor& v) const
  {
    v.visit(static_cast<const Array&>(*this));
  }

  virtual ~Array()
  {
    for(std::vector<Value*>::iterator it = this->begin();
        it != this->end(); ++it)
    {
      delete *it;
    }
  }
};

struct String : std::string, Value
{
  virtual void accept(Value::Visitor& v) const
  {
    v.visit(static_cast<const String&>(*this));
  }

  virtual ~String() { }
};

struct Number : Value
{
  double val;

  Number(std::string s)
  {
    std::stringstream ss(s);
    ss >> val;
  }

  virtual void accept(Value::Visitor& v) const
  {
    v.visit(static_cast<const Number&>(*this));
  }

  virtual ~Number() { }
};

struct Bool : Value
{
  bool val;

  Bool(bool b) : val(b) { }

  virtual void accept(Value::Visitor& v) const
  {
    v.visit(static_cast<const Bool&>(*this));
  }

  virtual ~Bool() { }
};

struct Null : Value
{
  virtual void accept(Value::Visitor& v) const
  {
    v.visit(static_cast<const Null&>(*this));
  }

  virtual ~Null() { }
};

/*std::ostream& operator<<(std::ostream& os, const Object& t)
{
  return os << t.val;
}*/

void print(std::ostream& out, Value& obj)
{
  struct V : Value::Visitor
  {
    V(std::ostream& o) : out(o), tab_level(0) { }

    void visit(const Object& obj)
    {
      out << "{" << std::endl;
      ++tab_level;
      for(Object::const_iterator it = obj.begin(); it != obj.end(); ++it)
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
    void visit(const Array& obj)
    {
      out << "[";
      for(Array::const_iterator it = obj.begin(); it != obj.end(); ++it)
      {
        (*it)->accept(*this);
        if(++it != obj.end())
          out << ",";
        --it;
      }
      out << "]";
    }
    void visit(const String& obj) { out << obj; }
    void visit(const Number& obj) { out << obj.val; }
    void visit(const Bool& obj) { out << (obj.val?"true":"false"); }
    void visit(const Null& obj) { out << "null"; }
    int tab_level;
    std::ostream& out;
  };
  V vis(out);
  obj.accept(vis);
}

Object* parse_object(std::istream&,JSON_checker&,int&);
Array* parse_array(std::istream&,JSON_checker&,int&);

Object parse(std::istream& in, int max_depth)
{
  JSON_checker jc = new_JSON_checker(max_depth);
  char cur;
  int line = 1;
  while(in >> cur)
  {
    if(JSON_checker_char(jc,cur))
    {
      if(cur == '\n')
        ++line;
      int state = jc->state;

      return *parse_object(in,jc,line);
    }else
    {
      std::cout << "Invalid JSON Syntax. at line: " << line
                << " character: " << cur << std::endl;
      exit(1);
    }
  }
}

Object* parse_object(std::istream& in, JSON_checker& jc, int& line)
{
  char cur;
  Object* document = new Object();
  std::string temp_key("");
  std::string temp_num("");
  Value * temp_value;
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
            Number* num = new Number(temp_num);
            temp_value = num;
            temp_num = "";
          }
          have_value = true;
          break;
        case 2: // Object
          temp_value = parse_object(in,jc,line);
          break;
        case 3: // Key
          if(temp_num.size() != 0)
          {
            Number* num = new Number(temp_num);
            temp_value = num;
            temp_num = "";
            have_value = true;
          }
          break;
        case 4: // Colon
          have_key = true;
          break;
        case 6: // Array
          temp_value = parse_array(in,jc,line);
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
          Bool* b_1 = new Bool(true);
          temp_value = b_1;
          break;
        }
        case 26: // False
        {
          Bool* b2 = new Bool(false);
          temp_value = b2;
          break;
        }
        case 29: // Null
          Null* n = new Null();
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

Array* parse_array(std::istream& in, JSON_checker& jc, int& line)
{
  char cur;
  Array* array = new Array();
  std::string temp_num("");
  std::string temp_string("");
  Value * temp_value;
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
            Number* num = new Number(temp_num);
            temp_value = num;
            temp_num = "";
          }
          have_value = true;
          break;
        case 2: // Object
          temp_value = parse_object(in,jc,line);
          break;
        case 5: // Value
          if(temp_num.size() != 0)
          {
            Number* num = new Number(temp_num);
            temp_value = num;
            temp_num = "";
            have_value = true;
          }
          break;
        case 6: // Array
          temp_value = parse_array(in,jc,line);
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
          Bool* b_1 = new Bool(true);
          temp_value = b_1;
          break;
        }
        case 26: // False
        {
          Bool* b2 = new Bool(false);
          temp_value = b2;
          break;
        }
        case 29: // Null
          Null* n = new Null();
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

} // namespace json

#endif // JSON_H
