#include <iostream>
#include "JSON_checker.h"
#include "json.h"

using namespace std;

void print_state(int s);
void print_mode(int s);
void test();

int main()
{
  test();

  /*
  JSON_checker jc = new_JSON_checker(20);
  char c;
  while(cin >> c)
  {
    if(!JSON_checker_char(jc,c))
    {
      cout << "syntax error" << endl;
      return 1;
    }
    cout << c << endl;
    print_state(jc->state);
  }
  if(!JSON_checker_done(jc))
  {
    cout << "syntax error" << endl;
    return 1;
  }
  //*/
}

void test()
{
  json::Object o = json::parse(cin,10);
  json::print(cout,o);
}

void print_mode(int s)
{
  switch(s)
  {
    case 0:
      cout << "MODE_ARRAY";
      break;
    case 1:
      cout << "MODE_DONE";
      break;
    case 2:
      cout << "MODE_KEY";
      break;
    case 3:
      cout << "MODE_OBJECT";
      break;
  }
}

void print_state(int s)
{
  switch(s)
  {
    case 0:
      cout << "GO:start" << endl;
      break;
    case 1:
      cout << "OK:ok" << endl;
      break;
    case 2:
      cout << "OB:object" << endl;
      break;
    case 3:
      cout << "KE:key" << endl;
      break;
    case 4:
      cout << "CO:colon" << endl;
      break;
    case 5:
      cout << "VA:value" << endl;
      break;
    case 6:
      cout << "AR:array" << endl;
      break;
    case 7:
      cout << "ST:string" << endl;
      break;
    case 8:
      cout << "ES:escape" << endl;
      break;
    case 9:
      cout << "U1:u1" << endl;
      break;
    case 10:
      cout << "U2:u2" << endl;
      break;
    case 11:
      cout << "U3:u3" << endl;
      break;
    case 12:
      cout << "U4:u4" << endl;
      break;
    case 13:
      cout << "MI:minus" << endl;
      break;
    case 14:
      cout << "ZE:zero" << endl;
      break;
    case 15:
      cout << "IN:integer" << endl;
      break;
    case 16:
      cout << "FR:fraction" << endl;
      break;
    case 17:
      cout << "E1:e" << endl;
      break;
    case 18:
      cout << "E2:ex" << endl;
      break;
    case 19:
      cout << "E3:exp" << endl;
      break;
    case 20:
      cout << "T1:tr" << endl;
      break;
    case 21:
      cout << "T2:tru" << endl;
      break;
    case 22:
      cout << "T3:true" << endl;
      break;
    case 23:
      cout << "F1:fa" << endl;
      break;
    case 24:
      cout << "F2:fal" << endl;
      break;
    case 25:
      cout << "F3:fals" << endl;
      break;
    case 26:
      cout << "F4:false" << endl;
      break;
    case 27:
      cout << "N1:nu" << endl;
      break;
    case 28:
      cout << "N2:nul" << endl;
      break;
    case 29:
      cout << "N3:null" << endl;
      break;
    case 30:
      cout << "NR_STATES" << endl;
      break;
  }
}
