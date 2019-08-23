#include "mcp3008Spi.h"
using namespace std;

//######################################################### Constructor
/** Default constructor. Set member variables to default values and 
then call spiOpen().
\author Jan */
MCP3008SPI::MCP3008SPI()
{
    this->mode = SPI_MODE_0 ;
    this->bitsPerWord = 8;
    this->speed = 8000000;
    this->spifd = -1;
 
    this->SpiOpen(std::string("/dev/spidev0.0"));
}
 
//######################################################### Overloaded Constructor
/** Overloaded constructor. Let user set member variables 
and then call spiOpen()
\author Jan */
MCP3008SPI::MCP3008SPI(std::string devspi, unsigned char spiMode, unsigned int spiSpeed, unsigned char spibitsPerWord)
{
    this->mode = spiMode ;
    this->bitsPerWord = spibitsPerWord;
    this->speed = spiSpeed;
    this->spifd = -1;
 
    this->SpiOpen(devspi);
    
    this->bufferPos = 0;
    for (int k=0;k<MCP3008_BUFFER_SIZE;k++)
    {
		this->buffer[k]=0;
	}
}
 
//######################################################### Destructor
/** Calls the spiClose() method.
\author Jan */
MCP3008SPI::~MCP3008SPI()
{
    this->SpiClose();
}
