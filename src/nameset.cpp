/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nameset.cpp,v 1.1 2001/11/06 16:18:31 bzfkocht Exp $"

#include <string.h>
#include "nameset.h"

namespace soplex
{

char NameSet_Name::deflt = '\0';

void NameSet::add(const char* str)
{
   Key key;
   add(key, str);
}

void NameSet::add(Key& key, const char* str)
{
   const NameSet_Name nstr (str);
   if (!hashtab.has(nstr))
   {
      if (size() + 1 > max())
      {
         assert(factor >= 1);
         reMax(int(factor*max() + 8));
      }

      if (memSize() + (int)strlen(str) >= memMax())
      {
         memPack();
         if (memSize() + (int)strlen(str) >= memMax())
         {
            assert(memFactor >= 1);
            memRemax(int(memFactor*memMax()) + 9 + strlen(str));
            assert(memSize() + (int)strlen(str) < memMax());
         }
      }

      char* tmp = &(mem[memused]);
      memused += strlen(str) + 1;
      strcpy(tmp, str);

      DataSet < NameSet_CharPtr > ::Key* keyptr = (DataSet < NameSet_CharPtr > ::Key*) & key;
      NameSet_CharPtr* name = set.create(*keyptr);
      name->name = tmp;
      list.append(name);
      NameSet_Name memname(name->name);
      hashtab.add(memname, key);
   }
}

void NameSet::add(const NameSet& set)
{
   for (int i = 0; i < set.num(); ++i)
   {
      NameSet_Name iname(set[i]);
      if (!hashtab.has(iname))
         add(set[i]);
   }
}

void NameSet::add(Key key[], const NameSet& set)
{
   for (int i = 0; i < set.num(); ++i)
   {
      NameSet_Name iname = set[i];
      if (!hashtab.has(iname))
         add(key[i], set[i]);
   }
}


void NameSet::remove(const char *str)
{
   const NameSet_Name nam(str);
   if (hashtab.has (nam))
   {
      DataSet < NameSet_CharPtr > ::Key* key = (DataSet < NameSet_CharPtr > ::Key*)hashtab.get (nam);
      hashtab.remove (nam);
      list.remove(&(set[*key]));
      set.remove(*key);
   }
}

void NameSet::remove(Key key)
{
   assert(has(key));
   DataSet < NameSet_CharPtr > ::Key* keyptr = (DataSet < NameSet_CharPtr > ::Key*) & key;
   const NameSet_Name nam = set[*keyptr].name;
   hashtab.remove (nam);
   list.remove(&(set[*keyptr]));
   set.remove(*keyptr);
}

void NameSet::remove(Key keys[], int n)
{
   for (int i = 0; i < n; ++i)
      remove(keys[i]);
}

void NameSet::remove(int nums[], int n)
{
   for (int i = 0; i < n; ++i)
      remove(nums[i]);
}

void NameSet::clear()
{
   set.clear();
   hashtab.clear();
   list.clear();
   memused = 0;
}



void NameSet::reMax(int newmax)
{
   hashtab.reMax (newmax);

   long delta = set.reMax(newmax);
   NameSet_CharPtr* first = list.first();

   if (delta != 0 && first != 0)
   {
      NameSet_CharPtr* last = (NameSet_CharPtr*)((char*)list.last() + delta);
      NameSet_CharPtr* name = (NameSet_CharPtr*)((char*)first + delta);
      first = name;
      for (; name != last; name = name->next())
         name->next() = (NameSet_CharPtr*)((char*)name->next() + delta);

      list = IsList < NameSet_CharPtr > (first, last);
   }
}

void NameSet::memRemax(int newmax)
{
   char* old = mem;
   long delta;

   memmax = (newmax < memSize()) ? memSize() : newmax;
   mem = (char*)realloc(mem, memmax * sizeof(char));
   if (mem == 0)
   {
      std::cerr << "ERROR: NameSet could not reallocate memory\n";
      exit(-1);
   }
   delta = mem - old;

   hashtab.clear ();

   for (NameSet_CharPtr* name = list.first(); name; name = list.next(name))
      name->name = &(name->name[delta]);

   for (int i = num() - 1; i >= 0; --i)
   {
      Key ikey = key(i);
      DataSet < NameSet_CharPtr > ::Key* keyptr = (DataSet < NameSet_CharPtr > ::Key*) & ikey;
      NameSet_Name nam(set[*keyptr].name);
      hashtab.add(nam, *keyptr);
   }
}

void NameSet::memPack()
{
   int i;

   hashtab.clear ();

   NameSet_CharPtr* name = list.first();
   for (memused = 0; name != NULL; name = list.next(name))
   {
      for (i = 0; (mem[memused + i] = name->name[i]) != 0; ++i)
        ;
      name->name = &(mem[memused]);
      memused += i + 1;
   }
   assert(memSize() <= memMax());

   for (i = num() - 1; i >= 0; --i)
   {
      Key ikey = key(i);
      DataSet < NameSet_CharPtr > ::Key* keyptr = (DataSet < NameSet_CharPtr > ::Key*) & ikey;
      NameSet_Name nam(set[*keyptr].name);
      hashtab.add(nam, *keyptr);
   }
}

static int hashFunction (const NameSet_Name* str)
{
   unsigned int res = 0;
   const char* sptr = str->name;

   while (*sptr)
   {
      res *= 65;
      res += *sptr++ - int('0');
      res %= 0x0fffffff;
   }
   return res;
}

NameSet& NameSet::operator=(const NameSet& rhs)
{
   if (max() < rhs.size())
      reMax(rhs.size());
   if (memMax() < rhs.memSize())
      memRemax(rhs.memSize());

   set = rhs.set;

   list.clear();
   hashtab.clear();
   for (int i = 0; i < set.num(); ++i)
   {
      list.append(&(set[i]));
      NameSet_Name iname(set[i].name);
      Key ikey = Key(set.key(i));
      hashtab.add(iname, ikey);
   }
   memPack();

   return *this;
}

NameSet::NameSet(const NameSet& org)
   : set(org.set)
      , hashtab(org.hashtab)
      , factor(org.factor)
      , memFactor(org.memFactor)
{
   memused = 0;
   memmax = org.memSize();
   mem = (char*)malloc(memmax * sizeof(char));
   if (mem == 0)
   {
      std::cerr << "ERROR: NameSet could not allocate memory\n";
      exit(-1);
   }

   list.clear();
   hashtab.clear();
   for (int i = 0; i < set.num(); ++i)
   {
      list.append(&(set[i]));
      NameSet_Name iname = set[i];
      Key key = Key(set.key(i));
      hashtab.add(iname, key);
   }
   memPack();
}

NameSet::NameSet(int max, int mmax, double fac, double memFac)
   : set(max)
      , hashtab(hashFunction, set.max(), 0, fac)
      , factor(fac)
      , memFactor(memFac)
{
   memused = 0;
   memmax = (mmax < 1) ? (8 * set.max() + 1) : mmax;
   mem = (char*)malloc(memmax * sizeof(char));
   if (mem == 0)
   {
      std::cerr << "ERROR: NameSet could not allocate memory\n";
      exit(-1);
   }
}

NameSet::~NameSet()
{
   free(mem);
}

int NameSet::isConsistent() const
{
   if (memused > memmax)
   {
      std::cerr << "Inconsistency detected in class NameSet!\n";
      return 0;
   }

   NameSet_CharPtr* next;
   for (NameSet_CharPtr *name = list.first(); name; name = next)
   {
      int len = strlen(name->name) + 1;
      if (&(name->name[len]) > &(mem[memused]))
      {
         std::cerr << "Inconsistency detected in class NameSet!\n";
         return 0;
      }

      next = list.next(name);
      if (next != 0 && &(name->name[len]) > next->name)
      {
         std::cerr << "Inconsistency detected in class NameSet!\n";
         return 0;
      }
   }

   return set.isConsistent() && list.isConsistent() && hashtab.isConsistent();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
