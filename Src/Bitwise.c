/*  Performs bit level operations on words for the Adrasteia program
    Copyright (C) 2006    Miles Hampson

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#include <stdio.h>
#include <stdlib.h>
//the wordlength (in bits) used to store unsigned ints
//(our bitword storage containers) on this machine
int WORDLENGTH;
//the actual wordlength of the strings in bit_strings
//this will be rounded up to the nearest 32 bit 
//word boundary for storage in bit_strings, the
//extra padding is discarded when bit slicing
int STRINGWORDLENGTH;
//the number of bits in each word to set to 1
int NUMSETBITS;
//An array of ints, each one representing a bit string
unsigned int *bit_strings;
//the number of segments each string has to be broken into
//due to it being larger than the wordlength
int bit_string_segments;
//the number of logical elements in the bit_string array
int bit_string_elms;
//the bits to set for the current code word
unsigned int *result;
//An array of horizontal slices of the bit_string array
unsigned int *bit_slices;
//the number of words needed to store a slice of all 
//records for a bit string array
int record_segments;
//An array of records that sucessfully match some input
unsigned int *dropset;
//the records that match each query piece
unsigned int *segments;
//number of sections
int sections;
//Array of all bit slice arrays
unsigned int **slices;
//the size of this array
int slice_size;

//Set up parameters for this database
void initialize(int word_len,int num_set_bits)
{
    WORDLENGTH = sizeof(unsigned int) * 8;
    STRINGWORDLENGTH = word_len;
    NUMSETBITS = num_set_bits;
    bit_string_segments = ((word_len-1)/WORDLENGTH)+1;
    if((result=malloc(sizeof(unsigned int)*bit_string_segments))==NULL)
    {
        printf("Insufficent memory for creation of result array\n");
        exit(1); 
    }
    bit_string_elms=0;
    record_segments=0;
    sections = 0;
    slice_size=0;
}

//Writes the completed bit_slices array to the specified location
void writeBitSlices(char *file_loc)
{
    int num_ints = record_segments*STRINGWORDLENGTH;
    FILE *fp;
    if((fp=fopen(file_loc, "wb"))==NULL) 
    {
        printf("Cannot open a slices file for writing.\n");
        exit(1);
    }
    if(fwrite(bit_slices, sizeof(unsigned int), num_ints, fp) != num_ints) 
    {
        printf("Error writing to a slices file.\n");
        exit(1);
    }
    fclose(fp);
    free(bit_slices);
    bit_slices=NULL;
}

//Load all slice files from disk to slice array
//Assumes that the idx files are consecutive between start and end
//and that they are all in the range 0 to 9
void loadBitSlices(int start_idx,int end_idx,int pos_to_replace,char *firstName)
{
    int i;
    int num_ints = record_segments*STRINGWORDLENGTH;
    FILE *fp;
    slice_size = (end_idx-start_idx)+1;
    if ((slices=malloc(sizeof(unsigned int *)*slice_size))==NULL)
    {
        printf("Not enough memory for slice array\n");
        exit(1); 
    }
    for(i=start_idx;i<=end_idx;i++)
    {
        firstName[pos_to_replace]=i + '0';
        if((fp=fopen(firstName, "rb"))==NULL) 
        {
            printf("Cannot open slice file for reading.\n");
            exit(1);
        }
        if ((slices[i-start_idx]=malloc(sizeof(unsigned int )*num_ints))==NULL)
        {
            printf("Not enough memory for slice array row\n");
            exit(1); 
        }
        if(fread(slices[i-start_idx],sizeof(unsigned int),num_ints,fp) != num_ints) 
        {
            if(feof(fp)) printf("Premature end of slice file (possibly due to out of date index files)\n");
            else printf("Slice file read error.");
            exit(1);
        }
        fclose(fp);
    }
}

//Resets indexing information about the database
//and returns the number of bit_string_elms
int resetDBInfo()
{
    int tmp_bse = bit_string_elms;
    record_segments=0;
    bit_string_elms=0;
    return tmp_bse;
}

//Sets indexing information for the database after processing
//Returns the previous number of record_segments
//if a 0 is given in the first parameter
int setDBInfo(int rec_segs,int bselms)
{
    if(rec_segs!=0)
        record_segments=rec_segs;
    bit_string_elms=bselms;
    return record_segments;
}

//Creates a new array of bit slices from the stored array of
//bit strings, and then frees the memory in this old array
void bitSlice()
{
    int j;
    int i;
    int bit_word_elms=((bit_string_elms-1)%WORDLENGTH)+1;
    if((bit_slices=realloc(bit_slices,sizeof(unsigned int)*
                            ++record_segments*STRINGWORDLENGTH))==NULL)
    {
        printf("bit_slices out of memory\n");
        exit(1); 
    } 
    for (i=0;i<STRINGWORDLENGTH;i++)
        bit_slices[(record_segments-1)*STRINGWORDLENGTH+i]=0;
    for(i=0;i<STRINGWORDLENGTH;i++)
        for(j=0;j<bit_word_elms;j++)
            if((1<<(i%WORDLENGTH)) & 
                    bit_strings[(j*bit_string_segments)+(i/WORDLENGTH)])
                bit_slices[((record_segments-1)*STRINGWORDLENGTH)+i] |= 1<<j;
    free(bit_strings);
}

//Allocates memory and stores information about a new record to be parsed
void allocateRecord()
{
    int i;
    int alloc_size;
    if(bit_string_elms%WORDLENGTH==0)
    {
        if(bit_string_elms!=0)
        {
            bitSlice();
        }
        alloc_size=bit_string_segments*WORDLENGTH;
        if ((bit_strings=malloc(sizeof(unsigned int)*alloc_size))==NULL)
        {
            printf("bit_strings out of memory\n");
            exit(1); 
        }
        for (i=0;i<alloc_size;i++)
            bit_strings[i]=0;
    }
    bit_string_elms+=1;
}

//Encode and store a piece of the database using
//the method of Superimposed Code Words
void encodeDBPiece(char *piece)
{
    int cur_bit;
    unsigned long hash = 5381;
    int i;
    for (i=0;i<bit_string_segments;i++)
        result[i]=0;
    while (i = *piece++)
        hash = ((hash << 5) + hash) + i;   
    srand(hash);
    for(i=0;i<NUMSETBITS;i++)
    {
         cur_bit = STRINGWORDLENGTH * (rand() / (RAND_MAX + 1.0));
         while((1<<(cur_bit%WORDLENGTH))&(result[cur_bit/WORDLENGTH]))
            cur_bit = STRINGWORDLENGTH * (rand() / (RAND_MAX + 1.0));
         result[cur_bit/WORDLENGTH] |= (1<<(cur_bit%WORDLENGTH));
         bit_strings[(((bit_string_elms-1)%WORDLENGTH)*bit_string_segments)
                            +cur_bit/WORDLENGTH] |= result[cur_bit/WORDLENGTH];
     }
}

//Encode and store a piece of a query using
//the method of Superimposed Code Words
void encodeQueryPiece(int seg_num,char *piece,int idx_num)
{
    int cur_bit;
    unsigned long hash = 5381;
    int i;
    int j;
    int seg_pos = seg_num*record_segments;
    while (i = *piece++)
        hash = ((hash << 5) + hash) + i;
    srand(hash);
    for (i=0;i<bit_string_segments;i++)
        result[i]=0;
    for(i=0;i<NUMSETBITS;i++)
    {
         cur_bit = STRINGWORDLENGTH * (rand() / (RAND_MAX + 1.0));
         while((1<<(cur_bit%WORDLENGTH))&(result[cur_bit/WORDLENGTH]))
            cur_bit = STRINGWORDLENGTH * (rand() / (RAND_MAX + 1.0));
         result[cur_bit/WORDLENGTH] |= 1 << cur_bit%WORDLENGTH; 
         for(j=0;j<record_segments;j++)
            segments[seg_pos+j] &= slices[idx_num][j*STRINGWORDLENGTH+cur_bit]; 
    }
}

void createSegments(int secs)
{
    int i;
    int section_size = record_segments*secs;
    sections = secs;
    if(segments!=NULL)free(segments);
    if((segments = malloc(sizeof(unsigned int)*section_size))==NULL) 
    {
        printf("Insufficent memory for creation of section array\n");
        exit(1); 
    } 
    for(i=0;i<section_size;i++)
        segments[i]=-1;
}

//Creates a set of the records that match when all segments of this query
//are or'ed with the previous dropset
void updateDropsetOr()
{
    int i;
    int j;
    int secs = record_segments*sections;
    if(dropset==NULL)
    {
        if((dropset = malloc(sizeof(unsigned int)*record_segments))==NULL)
        {
            printf("Insufficent memory for creation of dropset array\n");
            exit(1);
        }
        for(i=0;i<record_segments;i++)
            dropset[i]=0;
    }
    for(i=0;i<record_segments;i++)
        for(j=0;j<secs;j+=record_segments)
            dropset[i] |= segments[i+j];
}

//Creates a set of the records that match when all segments of this query
//are and'ed with the previous dropset
void updateDropsetAnd()
{
    int i;
    int j;
    int secs = record_segments*sections;
    unsigned int tmp;
    if(dropset==NULL)
    {
        if((dropset = malloc(sizeof(unsigned int)*record_segments))==NULL)
        {
            printf("Insufficent memory for creation of dropset array\n");
            exit(1);
        }
        for(i=0;i<record_segments;i++)
            dropset[i]=-1;
    }
    for(i=0;i<record_segments;i++)
    {
        tmp = 0;
        for(j=0;j<secs;j+=record_segments)
            tmp |= segments[i+j];
        dropset[i] &= tmp;
    }
}

//returns the position of the next drop, starting from the position
//given in pos (counting from 1). Returns 0 if reached end of recs
int nextDrop(int pos)
{
    int i;
    for(i=pos-1;i<bit_string_elms;i++)
        if((1<<(i%WORDLENGTH))&dropset[i/WORDLENGTH])return i+1;
    return 0;
}

void clearMatchInfo()
{
    if(dropset!=NULL)free(dropset);
    dropset=NULL;
    if(segments!=NULL)free(segments);
    segments=NULL;
}

void clearDatabase()
{
    int i;
    if(slices!=NULL)
    {
        for(i=0;i<slice_size;i++)
            free(slices[i]);
        free(slices);
    }
    slices=NULL;
    if(result!=NULL)free(result);
    result=NULL;
}

//Testing fucntion displaying an int as a binary string
void print_int_as_binary_string(unsigned int number)
{
    char *binary_string = malloc(WORDLENGTH+1);
    int i;
    binary_string[WORDLENGTH] = '\0';
    for(i=0;i<WORDLENGTH;i++)
    {
        binary_string[i]=((number&1)+'0');
        number=number>>1;
    }
    printf(binary_string);
}
 

