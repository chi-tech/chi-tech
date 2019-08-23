#include "chi_pi3.h"
#include <iostream>

//######################################################### Default Constructor
/**Constructor
\author Jan*/
CHI_PI3::CHI_PI3()
{
	this->spi = NULL;
	for (int k=0;k<=40;k++)
	{
		this->pinNumber[k]	= std::to_string(k);
		this->pinExported[k]	= false;
		//printf("Pin number: %s\n",pinNumber[k].c_str());
	}
	for (int k=0; k<8; k++)
	{
		this->bufferMCP3008SPI[k]=false;
	}
	
	uart0_filestream = -1;
}

//######################################################### Default Destructor
/**Constructor
\author Jan*/
CHI_PI3::~CHI_PI3()
{
	for (int k=0;k<=40;k++)
	{
		if (this->pinExported[k])
		{
			this->UnExportPin(k);
		}
	}
	
}
