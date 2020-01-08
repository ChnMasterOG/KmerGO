// coding = utf-8
// author: QiChen
// version: v3.0
// modification date: 2019/10/25

#include <iostream>
#include <stdlib.h>

#define MORE	        0	    //more than
#define LESS	        1	    //less than
#define EQUAL	        2	    //equal to
#define ERROR	        -1	    //error

//KLEN = K-mer length
//CLEN = capacity of loser tree
//a winner is in ls[0], losers are in others
int *ls, *CLEN, *KLEN;

extern "C"{

void Setup_parameters(int klen, int clen);
void Free_memory();
void Loser_Tree_Adjust(int s, char *buf);
void Loser_Tree_Build(char *buf);
int Loser_Tree_GetMin(void);

}

/* Setup parameters */
void Setup_parameters(int klen, int clen)
{
    CLEN = (int* )malloc(sizeof(int));
    KLEN = (int* )malloc(sizeof(int));
    ls =  (int* )malloc(sizeof(int) * (clen + 1));
    *CLEN = clen;
    *KLEN = klen;
}

/* Free the memory */
void Free_memory()
{
    free(CLEN);
    free(KLEN);
    free(ls);
}

/* Compare K-mer string code and return MORE, LESS or EQUAL */
char KmerCodeCompare(const char* p1, const char* p2)
{
	for (int i = 0; i < *KLEN; i++)
	{
		if (p1[i] > p2[i])	return MORE;
		if (p1[i] < p2[i])	return LESS;
		if (i == *KLEN - 1)	return EQUAL;
	}
	return ERROR;
}

/* Swap two variables */
void Swap(int* a, int* b)
{
	int tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

/* Adjust the loser tree and s is the index */
void Loser_Tree_Adjust(int s, char *buf)
{
	int t = (s + *CLEN) / 2;		//get the father node of s
	while (t > 0) {					//compare with the parent node until the root node of the loser tree
		if (KmerCodeCompare(buf + s * (*KLEN + 1), buf + ls[t] * (*KLEN + 1)) == MORE)
		{
			Swap(&s, &ls[t]);		//swap s and ls[t]
		}
		t /= 2;
	}
	ls[0] = s;						//the last winner is in ls[0]
}

/* Build a loser tree */
void Loser_Tree_Build(char *buf)
{
	int i;
	for (i = 0; i < *CLEN + 1; ++i)
		ls[i] = *CLEN;
	for (i = 0; i < *CLEN; ++i)
		Loser_Tree_Adjust(i, buf);	//adjust the position one by one
}

/* Get the minimum position*/
int Loser_Tree_GetMin(void)
{
	return ls[0];
}
