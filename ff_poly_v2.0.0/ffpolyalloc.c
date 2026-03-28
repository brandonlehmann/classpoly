/*
    Copyright 2017 Andrew V. Sutherland

    This file is part of ff_poly.

    ff_poly is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2 of the License.

    ff_poly is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ff_poly.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include "ffpolyalloc.h"

struct ff_poly_stack_context *_ff_sctx;

void ff_poly_stack_clear (void)
{
    if ( ! _ff_sctx ) return;
    if ( _ff_sctx->next_entry ) fprintf (stderr, "WARNING: ff_poly_stack_clear called with %d entries still allocated (%d offstack)\n", _ff_sctx->next_entry, _ff_sctx->offstack_entries);
    if ( _ff_sctx->offstack_entries ) {
        for ( int i = 0 ; i < _ff_sctx->next_entry ; i++ ) if ( ! _ff_sctx->tab[i].onstack ) mem_free (_ff_sctx->tab[i].ptr);
    }
    free (_ff_sctx->base);
    free (_ff_sctx->tab);
    free (_ff_sctx);
    _ff_sctx = 0;
}

void ff_poly_stack_init (int size, int entries)
{
    if ( size < FF_POLY_DEFAULT_STACK_SIZE ) size = FF_POLY_DEFAULT_STACK_SIZE;
    if ( entries < FF_POLY_DEFAULT_STACK_ENTRIES ) entries = FF_POLY_DEFAULT_STACK_ENTRIES;
    // Avoid re-allocating stack if we can
    if ( _ff_sctx && (_ff_sctx->end-_ff_sctx->base) >= size && _ff_sctx->entries >= entries && ! _ff_sctx->offstack_entries) { _ff_sctx->next_entry = 0; return; }
    ff_poly_stack_clear ();
    _ff_sctx = malloc (sizeof(*_ff_sctx));
    _ff_sctx->base = _ff_sctx->next = malloc (size*sizeof(ff_t));
    _ff_sctx->end = _ff_sctx->base + size;
    _ff_sctx->tab = malloc (entries*sizeof(*_ff_sctx->tab));
    _ff_sctx->entries = entries;
    _ff_sctx->next_entry = 0;
    _ff_sctx->offstack_entries = 0;
}
