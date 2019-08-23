#include "mcp3008Spi.h"
using namespace std;

 
//######################################################### Exchange with SPI
/**This function uses a fundamental data structure to get the sample
from the chip.
\author Jan */
int MCP3008SPI::SpiWriteRead( unsigned char *data, int length)
{
  	struct spi_ioc_transfer spi[length];
  	int i = 0;
 	int retVal = -1; 
  	bzero(spi, sizeof spi); // ioctl struct must be zeroed 
 
	// one spi transfer for each byte
 
  	for (i = 0 ; i < length ; i++)
	{
 
    		spi[i].tx_buf        = (unsigned long)(data + i); // transmit from "data"
    		spi[i].rx_buf        = (unsigned long)(data + i) ; // receive into "data"
    		spi[i].len           = sizeof(*(data + i)) ;
    		spi[i].delay_usecs   = 0 ;
    		spi[i].speed_hz      = this->speed ;
    		spi[i].bits_per_word = this->bitsPerWord ;
    		spi[i].cs_change = 0;
	}
 
 	retVal = ioctl (this->spifd, SPI_IOC_MESSAGE(length), &spi) ;
 
 	if(retVal < 0)
	{
    		perror("Problem transmitting spi data..ioctl");
    		exit(1);
	}
 
	return retVal;
}


//######################################################### Read a Channel
/** Reads the chip's indicated channel.
\author Jan */
int MCP3008SPI::ReadChannel(int channel)
{
	unsigned char data[3];
	int a2dVal = 0;

	data[0] = 1;  //  first byte transmitted -> start bit
        data[1] = 0b10000000 |( ((channel & 7) << 4)); // second byte transmitted -> (SGL/DIF = 1, D2=D1=D0=0)
        data[2] = 0; // third byte transmitted....don't care

        this->SpiWriteRead(data, sizeof(data) );

        a2dVal = 0;
        a2dVal = (data[1]<< 8) & 0b1100000000; //merge data[1] & data[2] to get result
        a2dVal |=  (data[2] & 0xff);

	return a2dVal;
}
 

