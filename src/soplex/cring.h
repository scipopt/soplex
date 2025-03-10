/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
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
