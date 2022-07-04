#ifndef COMPRESSION_H
#define COMPRESSION_H

unsigned int code_compression(unsigned int, int, int, int *);/*declare code_compression*/

unsigned int CR_format(unsigned int, unsigned int, unsigned int, unsigned int); /*declare CR_format*/
unsigned int CI_format(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int); /*declare CI_format*/
unsigned int CL_format(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int); /*declare CL_format*/
unsigned int CS_format_1(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int); /*declare CS_format*/
unsigned int CS_format_2(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int); 
unsigned int CB_format_1(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int); /*declare CB_format*/
unsigned int CB_format_2(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int); 
unsigned int CJ_format(unsigned int, unsigned int, unsigned int); /*declare CJ_format*/
#endif