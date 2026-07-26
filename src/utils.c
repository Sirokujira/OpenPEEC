/*
utils.c

utilities (tokenize は OpenFDTD sol/utils.c から移植)
*/

#include "peec.h"

// tokenize a string
int tokenize(char *str, const char *tokensep, char *token[], int maxtoken)
{
	if ((str == NULL) || !maxtoken) return 0;

	char *thistoken = strtok(str, tokensep);

	int   count;
	for (count = 0; (count < maxtoken) && (thistoken != NULL); ) {
		token[count++] = thistoken;
		thistoken = strtok(NULL, tokensep);
	}

	token[count] = NULL;

	return count;
}
