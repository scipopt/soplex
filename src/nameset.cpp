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
#pragma ident "@(#) $Id: nameset.cpp,v 1.15 2002/01/22 14:17:15 bzfkocht Exp $"

#include <string.h>
#include "real.h"
#include "nameset.h"
#include "spxalloc.h"

namespace soplex
{
const char NameSet::Name::deflt = '\0';

void NameSet::add(const char* str)
{
   DataKey k;
   add(k, str);
}

void NameSet::add(DataKey& p_key, const char* str)
{
   const Name nstr (str);

   if (!hashtab.has(nstr))
   {
      if (size() + 1 > max())
      {
         assert(factor >= 1);
         reMax(int(factor*max() + 8));
      }

      if (memSize() + int(strlen(str)) >= memMax())
      {
         memPack();
         if (memSize() + int(strlen(str)) >= memMax())
         {
            assert(memFactor >= 1);
            memRemax(int(memFactor*memMax()) + 9 + int(strlen(str)));
            assert(memSize() + int(strlen(str)) < memMax());
         }
      }

      char* tmp = &(mem[memused]);
      memused += int(strlen(str)) + 1;
      strcpy(tmp, str);

      CharPtr* name = set.create(p_key);
      name->name = tmp;
      list.append(name);
      Name memname(name->name);
      hashtab.add(memname, p_key);
   }
}

void NameSet::add(const NameSet& p_set)
{
   for (int i = 0; i < p_set.num(); ++i)
   {
      Name iname(p_set[i]);
      if (!hashtab.has(iname))
         add(p_set[i]);
   }
}

void NameSet::add(DataKey p_key[], const NameSet& p_set)
{
   for (int i = 0; i < p_set.num(); ++i)
   {
      Name iname = Name(p_set[i]);
      if (!hashtab.has(iname))
         add(p_key[i], p_set[i]);
   }
}


void NameSet::remove(const char *str)
{
   const Name nam(str);
   if (hashtab.has (nam))
   {
      DataKey* hkey = hashtab.get(nam);
      hashtab.remove (nam);
      list.remove(&(set[*hkey]));
      set.remove(*hkey);
   }
}

void NameSet::remove(DataKey p_key)
{
   assert(has(p_key));
   const Name nam = set[p_key].name;
   hashtab.remove (nam);
   list.remove(&(set[p_key]));
   set.remove(p_key);
}

void NameSet::remove(DataKey keys[], int n)
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

   ptrdiff_t delta = set.reMax(newmax);
   CharPtr* first = list.first();

   if (delta != 0 && first != 0)
   {
      CharPtr* last = ( reinterpret_cast<CharPtr*>
                                (reinterpret_cast<char*>(list.last()) + delta) );
      CharPtr* name = ( reinterpret_cast<CharPtr*>
                                (reinterpret_cast<char*>(first) + delta) );
      first = name;
      for (; name != last; name = name->next())
         name->next() = reinterpret_cast<CharPtr*>
            (reinterpret_cast<char*>(name->next()) + delta);

      list = IsList < CharPtr > (first, last);
   }
}

void NameSet::memRemax(int newmax)
{
   char* old = mem;
   ptrdiff_t delta;

   memmax = (newmax < memSize()) ? memSize() : newmax;
   spx_realloc(mem, memmax);

   delta = mem - old;

   hashtab.clear ();

   /* update pointers to new targets */
   for (CharPtr* name = list.first(); name; name = list.next(name))
      name->name += delta;

   for (int i = num() - 1; i >= 0; --i)
   {
      Name nam(set[key(i)].name);
      hashtab.add(nam, key(i));
   }
}

void NameSet::memPack()
{
   int i;

   hashtab.clear ();

   CharPtr* name = list.first();
   for (memused = 0; name != 0; name = list.next(name))
   {
      for (i = 0; (mem[memused + i] = name->name[i]) != 0; ++i)
        ;
      name->name = &(mem[memused]);
      memused += i + 1;
   }
   assert(memSize() <= memMax());

   for (i = num() - 1; i >= 0; --i)
   {
      Name nam(set[key(i)].name);
      hashtab.add(nam, key(i));
   }
}

/// returns the hash value of the name.
int NameSetNameHashFunction(const NameSet::Name* str)
{
   unsigned int res = 0;
   const char* sptr = str->name;

   while(*sptr != '\0')
   {
      res *= 65;
      res += *sptr++ - int('0');
      res %= 0x0fffffff;
   }
   return res;
}

NameSet& NameSet::operator=(const NameSet& rhs)
{
   if (this != &rhs)
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
         Name iname(set[i].name);
         DataKey ikey = DataKey(set.key(i));
         hashtab.add(iname, ikey);
      }
      memPack();
   }
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
   spx_alloc(mem, memmax);

   list.clear();
   hashtab.clear();
   for (int i = 0; i < set.num(); ++i)
   {
      list.append(&(set[i]));
      Name iname = set[i];
      DataKey k = DataKey(set.key(i));
      hashtab.add(iname, k);
   }
   memPack();
}

NameSet::NameSet(int p_max, int mmax, Real fac, Real memFac)
   : set(p_max)
   , hashtab(NameSetNameHashFunction, set.max(), 0, fac)
   , factor(fac)
   , memFactor(memFac)
{
   memused = 0;
   memmax = (mmax < 1) ? (8 * set.max() + 1) : mmax;
   spx_alloc(mem, memmax);
}

NameSet::~NameSet()
{
   spx_free(mem);
}

bool NameSet::isConsistent() const
{
   if (memused > memmax)
      return MSGinconsistent("NameSet");

   CharPtr* next;

   for (CharPtr *name = list.first(); name; name = next)
   {
      int len = int(strlen(name->name)) + 1;

      if (&(name->name[len]) > &(mem[memused]))
         return MSGinconsistent("NameSet");

      next = list.next(name);

      if (next != 0 && &(name->name[len]) > next->name)
         return MSGinconsistent("NameSet");
   }
   return set.isConsistent() && list.isConsistent() && hashtab.isConsistent();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
