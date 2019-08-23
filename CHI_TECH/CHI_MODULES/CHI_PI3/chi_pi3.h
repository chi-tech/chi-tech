#ifndef CHI_PI3_H
#define CHI_PI3_H
#include<string>
#include"mcp3008spi/mcp3008Spi.h"

#define PI3_READ   0
#define PI3_WRITE  1

#define PI3_LOW  0
#define PI3_HIGH  1

#define PI3_BAUD_1200  0
#define PI3_BAUD_2400  1
#define PI3_BAUD_4800  2
#define PI3_BAUD_9600  3
#define PI3_BAUD_19200 4
#define PI3_BAUD_38400 5


//######################################################### CLASS DEF
/**Object for controlling Raspberry Pi3 GPIO.

\author Jan*/
class CHI_PI3
{
//=========================== Attributes
public:
	std::string pinNumber[41];
	MCP3008SPI* spi;
	bool        bufferMCP3008SPI[8];
private:
	bool		pinExported[41];
	int         uart0_filestream;
//=========================== Methods
public:
	//00 ConstrDestr
			CHI_PI3();
			~CHI_PI3();
	//01 Export Control
	bool		ExportPin(int pinNum);
	bool		UnExportPin(int pinNum);
	bool		SetPinMode(int pinNum, int mode);
	bool		SetPinValue(int pinNum, int value);
	int		GetPinValue(int pinNum);
	//02 Serial Communication
	bool        InitializeSerial(int baudrate);
	bool        SerialWrite(char* message);
	bool        SerialRead(char* message);
};

#endif
