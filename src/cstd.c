/*
    Copyright (c) 2007-2014 Andrew V. Sutherland

    This file is part of classpoly.

    classpoly is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2 of the License.

    classpoly is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with classpoly.  If not, see <http://www.gnu.org/licenses/>.
*/

int dbg_level;
unsigned long _cstd_seed;

int ui_qsort_cmp (const void *a, const void *b)
{
    if ( *((unsigned long *)a) < *((unsigned long*)b) ) return -1;
    if ( *((unsigned long *)a) > *((unsigned long*)b) ) return 1;
    return 0;
}

int dbl_qsort_cmp (const void *a, const void *b)
{
    if ( *((double *)a) < *((double*)b) ) return -1;
    if ( *((double *)a) > *((double*)b) ) return 1;
    return 0;
}
