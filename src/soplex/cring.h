/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _CRING_H_
#define _CRING_H_

/***************************************************************
                    Double linked ring
 ***************************************************************/

#define initDR(ring)    ((ring).prev = (ring).next = &(ring))

#define init2DR(elem, ring)                                     \
{                                                               \
(elem).next = (ring).next;                                 \
(elem).next->prev = &(elem);                               \
(elem).prev = &(ring);                                     \
(ring).next = &(elem);                                     \
}

#define removeDR(ring)                                          \
{                                                               \
(ring).next->prev = (ring).prev;                           \
(ring).prev->next = (ring).next;                           \
}

#endif // _CRING_H_
