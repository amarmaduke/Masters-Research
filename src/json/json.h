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

	Number() : val(0) { }

  Number(std::string s)
  {
    std::stringstream ss(s);
    ss >> val;
  }

	Number(double d) : val(d) { }

  virtual void accept(Value::Visitor& v) const
  {
    v.visit(static_cast<const Number&>(*this));
  }

  virtual ~Number() { }
};

struct Bool : Value
{
  bool val;

	Bool() : val(false) { }

  Bool(bool b) : val(b) { }

  virtual void accept(Value::Visitor& v) const
  {
    v.visit(static_cast<const Bool&>(*this));
  }

  virtual ~Bool() { }
};

struct Null : Value
{
	Null() { }

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

void print(std::ostream&,Value&);

Object parse(std::istream&,int);

Object* parse_object(std::istream&,JSON_checker&,int&);

Array* parse_array(std::istream&,JSON_checker&,int&);

} // namespace json

#endif // JSON_H
