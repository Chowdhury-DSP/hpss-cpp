#pragma once

#include <hpss/util/memory_arena.hpp>

namespace hpss
{
// adapted from: https://ideone.com/8VVEa (MIT license)

typedef float Item;
struct Mediator
{
    Item* data; //circular queue of values
    int* pos; //index into `heap` for each value
    int* heap; //max/median/min heap holding indexes into `data`.
    int N; //allocated size.
    int idx; //position in circular queue
    int ct; //count of items in queue
};

int MediatorSizeBytes (int nItems);

//creates new Mediator: to calculate `nItems` running median.
Mediator* MediatorNew (Memory_Arena& arena, int nItems);

//Inserts item, maintains median in O(lg nItems)
void MediatorInsert (Mediator* m, Item v);

//returns median item (or average of 2 when item count is even)
Item MediatorMedian (Mediator* m);
} // namespace hpss
