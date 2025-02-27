// ParallaxBA.cpp : Defines the exported functions for the DLL application.

#include "stdafx.h"
#include "SBA.h"
#include "../../../src/ParallaxBA/SBAImp.h"

SBAapi ISBA*	newSBA()
{
	return (ISBA*)(new CSBA );
}

SBAapi void    freeSBA( ISBA* ptr )
{
	free( ptr );
}