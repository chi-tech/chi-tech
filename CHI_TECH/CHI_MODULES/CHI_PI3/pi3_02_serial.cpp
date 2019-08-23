#include "chi_pi3.h"
#include <iostream>
#include <unistd.h>			//Used for UART
#include <fcntl.h>			//Used for UART
#include <termios.h>		//Used for UART

//#################################################################### Initialize Serial
/** \brief Initializes the Pi's UART serial communication.

\param baudrate   int Can be any of the <I>Serial Options</I>.

### Serial options
PI3_BAUD_1200  0
PI3_BAUD_2400  1
PI3_BAUD_4800  2
PI3_BAUD_9600  3
PI3_BAUD_19200 4
PI3_BAUD_38400 5

\example
CHI_PI3 chipie;
chipie.InitializeSerial(PI3_BAUD_9600);
\endexample

\author Jan
*/
bool CHI_PI3::InitializeSerial(int baudrate)
{
	
	uart0_filestream = open("/dev/serial0", O_RDWR | O_NOCTTY );		//Open in non blocking read/write mode
	if (uart0_filestream == -1)
	{
		printf("Error - Unable to open UART.  Ensure it is not in use by another application\n");
		return false;
	}
	
	//CONFIGURE THE UART
	//The flags (defined in /usr/include/termios.h - see http://pubs.opengroup.org/onlinepubs/007908799/xsh/termios.h.html):
	//	Baud rate:- B1200, B2400, B4800, B9600, B19200, B38400, B57600, B115200, B230400, B460800, B500000, B576000, B921600, B1000000, B1152000, B1500000, B2000000, B2500000, B3000000, B3500000, B4000000
	//	CSIZE:- CS5, CS6, CS7, CS8
	//	CLOCAL - Ignore modem status lines
	//	CREAD - Enable receiver
	//	IGNPAR = Ignore characters with parity errors
	//	ICRNL - Map CR to NL on input (Use for ASCII comms where you want to auto correct end of line characters - don't use for bianry comms!)
	//	PARENB - Parity enable
	//	PARODD - Odd parity (else even)
	struct termios options;
	speed_t inputBaud = B9600;
	
	switch (baudrate)
	{
		case PI3_BAUD_1200:  inputBaud = B1200;
		case PI3_BAUD_2400:  inputBaud = B2400;
		case PI3_BAUD_4800:  inputBaud = B4800;
		case PI3_BAUD_9600:  inputBaud = B9600;
		case PI3_BAUD_19200: inputBaud = B19200;
		case PI3_BAUD_38400: inputBaud = B38400;
		default:             inputBaud = B9600;
	}
	
	
	tcgetattr(uart0_filestream, &options);
	
	options.c_cflag = inputBaud | CS8 | CLOCAL | CREAD;		//<Set baud rate
	options.c_iflag = IGNPAR;
	options.c_oflag = 0;
	options.c_lflag = 0;
	
	tcflush(uart0_filestream, TCIFLUSH);
	tcsetattr(uart0_filestream, TCSANOW, &options);
	
	return true;
}

//#################################################################### Write to serial
/** \brief Writes a message to the serial port.
\param message char Message to be sent.*/
bool CHI_PI3::SerialWrite(char* message)
{
	char buf[]="hello world";
char buf2[11];

	
	if (uart0_filestream != -1)
	{
		int count = write(uart0_filestream, buf, 11);
		count=read(uart0_filestream,buf2,11);buf2[11]=0;
		//int count = write(uart0_filestream, message, 3);		//Filestream, bytes to write, number of bytes to write
		if (count <= 0)
		{
			printf("UART TX error\n");
			return false;
		}
	}
	return true;
}

//#################################################################### Read from serial
/** \brief Reads a message from the serial port.
\param message char Message read.*/
bool CHI_PI3::SerialRead(char* message)
{
	//----- CHECK FOR ANY RX BYTES -----
	if (uart0_filestream != -1)
	{
		// Read up to 255 characters from the port if they are there
		unsigned char rx_buffer[256];
		int rx_length = read(uart0_filestream, (void*)rx_buffer, 255);		//Filestream, buffer to store in, number of bytes to read (max)
		if (rx_length < 0)
		{
			return false;
		}
		else if (rx_length == 0)
		{
			//No data waiting
			return false;
		}
		else
		{
			//Bytes received
			rx_buffer[rx_length] = '\0';
			printf("%i bytes read : %s\n", rx_length, rx_buffer);
		}
	}
	return true;
}
