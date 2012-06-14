%module Bitwise
%{
extern void allocateRecord();
extern void initialize(int word_len,int num_set_bits);
extern int writeBitSlices(char *file_loc);
extern void loadBitSlices(int start_idx,int end_idx,int pos_to_replace,char *firstName);
extern int resetDBInfo();
extern int setDBInfo(int rec_segs,int bselms);
extern void encodeDBPiece(char *piece);
extern void bitSlice();
extern void encodeQueryPiece(int seg_num,char *piece,int slice_file_pos);
extern void createSegments(int sections);
extern void updateDropsetOr();
extern void updateDropsetAnd();
extern int nextDrop(int pos);
extern void clearMatchInfo();
extern void clearDatabase();
%}
extern void allocateRecord();
extern void initialize(int word_len,int num_set_bits);
extern int writeBitSlices(char *file_loc);
extern void loadBitSlices(int start_idx,int end_idx,int pos_to_replace,char *firstName);
extern int resetDBInfo();
extern int setDBInfo(int rec_segs,int bselms);
extern void encodeDBPiece(char *piece);
extern void bitSlice();
extern void encodeQueryPiece(int seg_num,char *piece,int slice_file_pos);
extern void createSegments(int sections);
extern void updateDropsetOr();
extern void updateDropsetAnd();
extern int nextDrop(int pos);
extern void clearMatchInfo();
extern void clearDatabase();
