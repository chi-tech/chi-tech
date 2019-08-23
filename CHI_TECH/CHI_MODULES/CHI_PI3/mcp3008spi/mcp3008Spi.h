/***********************************************************************
 * This header file contains the mcp3008Spi class definition.
 * Its main purpose is to communicate with the MCP3008 chip using
 * the userspace spidev facility.
 * The class contains four variables:
 * mode        -> defines the SPI mode used. In our case it is SPI_MODE_0.
 * bitsPerWord -> defines the bit width of the data transmitted.
 *        This is normally 8. Experimentation with other values
 *        didn't work for me
 * speed       -> Bus speed or SPI clock frequency. According to
 *                https://projects.drogon.net/understanding-spi-on-the-raspberry-pi/
 *            It can be only 0.5, 1, 2, 4, 8, 16, 32 MHz.
 *                Will use 1MHz for now and test it further.
 * spifd       -> file descriptor for the SPI device
 *
 * The class contains two constructors that initialize the above
 * variables and then open the appropriate spidev device using spiOpen().
 * The class contains one destructor that automatically closes the spidev
 * device when object is destroyed by calling spiClose().
 * The spiWriteRead() function sends the data "data" of length "length"
 * to the spidevice and at the same time receives data of the same length.
 * Resulting data is stored in the "data" variable after the function call.
 * ****************************************************************************/
#ifndef MCP3008SPI_H
#define MCP3008SPI_H
     
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/spi/spidev.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string>
#include <iostream>
 
#define MCP3008_BUFFER_SIZE 1000
 
class MCP3008SPI
{
	public:
		int             buffer[MCP3008_BUFFER_SIZE];
    	int             bufferPos;
private:
    	unsigned char 	mode;
    	unsigned char 	bitsPerWord;
    	unsigned int 	speed;
    	int 			spifd;
    	
    	
     
public:
	//00 ConstrDestr
    		MCP3008SPI();
    		MCP3008SPI(std::string devspi, unsigned char spiMode, unsigned int spiSpeed, unsigned char spibitsPerWord);
    		~MCP3008SPI();

	//01 OpenClose
private:
	int 	SpiOpen(std::string devspi);
    	int 	SpiClose();
	
	//02 Read
    	int 	SpiWriteRead( unsigned char *data, int length);
public:
    	int 	ReadChannel(int channel);
};
 
#endif
