#include "circular_buffer.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void ringbuffer_clear(ringbuffer_t *buffer, uint32_t size)
{
    buffer->S=size;

    uint32_t q = 0;
    for ( q = 0; q < size; q++)
    {
        buffer->m_buffer[q] = 0.0f;
    }

    buffer->m_size = 0;
    buffer->m_front = 0;
    buffer->m_back  = buffer->S - 1;
    buffer->power = 0.0f;
}

void ringbuffer_push(ringbuffer_t *buffer)
{
    buffer->m_back = (buffer->m_back + 1) % buffer->S;

    if(buffer->m_size == buffer->S)
    {
        buffer->m_front = (buffer->m_front + 1) % buffer->S;
    }
    else
    {
        buffer->m_size++;
    }
}
    
void ringbuffer_push_sample(ringbuffer_t *buffer, const float x)
{
    ringbuffer_push(buffer);
    buffer->m_buffer[buffer->m_back] = x;
}
    
void ringbuffer_pop(ringbuffer_t *buffer)
{
    if(buffer->m_size > 0 ) 
    {
        buffer->m_size--;
        buffer->m_front = (buffer->m_front + 1) % buffer->S;
    }
}
    
void ringbuffer_back_erase(ringbuffer_t *buffer, const uint32_t n)
{
    if(n >= buffer->m_size)
    {
        ringbuffer_clear(buffer, buffer->S);
    }
    else 
    {
        buffer->m_size -= n;
        buffer->m_back = (buffer->m_front + buffer->m_size - 1) % buffer->S;
    }
}
    
void ringbuffer_front_erase(ringbuffer_t *buffer, const uint32_t n)
{
    if(n >= buffer->m_size)
    {
        ringbuffer_clear(buffer, buffer->S);
    }
    else 
    {
        buffer->m_size -= n;
        buffer->m_front = (buffer->m_front + n) % buffer->S;
    }
}

int ringbuffer_peek_index(ringbuffer_t *buffer)
{
	uint32_t peek_index = 0;
	float peek_value = 0;
	
	uint32_t i = 0;
	for (i = 0; i < buffer->S; i++)
	{
		if (peek_value < buffer->m_buffer[i])
		{
			peek_value = buffer->m_buffer[i];
			peek_index = i;
		}
	}

	return peek_index;
}

float ringbuffer_push_and_calculate_power(ringbuffer_t *buffer, const float input)
{
    const float pow = (float)sqrt(input * input) * (1.0f / buffer->S);

    if (buffer->m_size < buffer->S)
    {
        //remove old sample and add new one to windowPower
        buffer->power += pow;
        ringbuffer_push_sample(buffer, pow);
    }
    else
    {
        //remove old sample and add new one to windowPower
        buffer->power += pow - ringbuffer_front(buffer);
        ringbuffer_pop(buffer);
        ringbuffer_push_sample(buffer, pow);
    }

    return buffer->power;
}

float ringbuffer_front(ringbuffer_t *buffer)
{ 
	return buffer->m_buffer[buffer->m_front]; 
}

float ringbuffer_back(ringbuffer_t *buffer)
{ 
	return buffer->m_buffer[buffer->m_back]; 
}

float ringbuffer_get_relative_val(ringbuffer_t *buffer, uint32_t index)
{ 
    if (buffer->m_back + index >= buffer->S)
        return buffer->m_buffer[(buffer->m_back + index) % (buffer->S)];
    else
	   return buffer->m_buffer[buffer->m_back + index]; 
}

float ringbuffer_get__direct_index_val(ringbuffer_t *buffer, uint32_t index)
{ 
    if (index >= buffer->S)
        return 0;

    return buffer->m_buffer[index]; 
}

int ringbuffer_empty(ringbuffer_t *buffer)
{ 
	return buffer->m_size == 0; 
}

int ringbuffer_full(ringbuffer_t *buffer)
{ 
	return buffer->m_size == buffer->S; 
}

float * ringbuffer_get_first_pointer(ringbuffer_t *buffer)
{
    return &buffer->m_buffer[buffer->m_front];
}
