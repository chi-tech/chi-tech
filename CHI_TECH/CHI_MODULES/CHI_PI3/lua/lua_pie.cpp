#include "../../../CHI_LUA/chi_lua.h"
#include "../chi_pi3.h"

extern CHI_PI3 chiPie;

/** \defgroup LuaPie Raspberry Pie
 * \ingroup LuaGeneralUtilities*/

//########################################################## Export Pin
/** \brief Exports a GPIO pin.

\param pinNumber int GPIO pin number to be exported.

\ingroup LuaPie
\author Jan*/
int chiPieExportPin(lua_State *L)
{
	int pinNum = lua_tonumber(L,1);
	
	chiPie.ExportPin(pinNum);
	
	return 0;
}

//########################################################## Set Pin Mode
/** \brief Sets a pin to READ or WRITE.

\param pinNumber int GPIO pin number.
\param mode int 0=Read, 1=Write.

\ingroup LuaPie
\author Jan*/
int chiPieSetPinMode(lua_State *L)
{
	int pinNum = lua_tonumber(L,1);
	int mode = lua_tonumber(L,2);
	
	chiPie.SetPinMode(pinNum,mode);
	
	return 0;
}

//########################################################## Set Pin Value
/** \brief Sets a pin to HIGH or LOW.

\param pinNumber int GPIO pin number.
\param mode int 0=Low, 1=High

\ingroup LuaPie
\author Jan*/
int chiPieSetPinValue(lua_State *L)
{
	int pinNum = lua_tonumber(L,1);
	int value = lua_tonumber(L,2);
	
	chiPie.SetPinValue(pinNum,value);
	
	return 0;
}

//########################################################## Get Pin Value
/** \brief Reads a pin.

\param pinNumber int GPIO pin number.

\ingroup LuaPie
\author Jan*/
int chiPieGetPinValue(lua_State *L)
{
	int pinNum = lua_tonumber(L,1);
	
	int value = 0;
	value = chiPie.GetPinValue(pinNum);
	
	lua_pushnumber(L,value);
	return 1;
}

//########################################################## Init SPI
/** \brief Initializes SPI.


\ingroup LuaPie
\author Jan*/
int chiPieInitSPI(lua_State *L)
{
	chiPie.spi = new MCP3008SPI;
	return 0;
}

//########################################################## Read SPI unbuffered
/** \brief Reads an SPI MCP3008 chip unbuffered.

\param channelNumber int Channel 0-7 to be read off the MCP3008 chip.

\ingroup LuaPie
\author Jan*/
int chiPieReadSPIChannel(lua_State *L)
{
	int channel = lua_tonumber(L,1);
	int value = 0;
	value = chiPie.spi->ReadChannel(channel);
	
	lua_pushnumber(L,value);
	return 1;
}

//########################################################## Set SPI buffer
/** \brief Sets an SPI MCP3008 chip to either be buffered or not.

\param channelNumber int Channel 0-7 to be read off the MCP3008 chip.
\param bufferFlag bool true=Is buffered, false= Not buffered.

\ingroup LuaPie
\author Jan*/
int chiPieSetSPIBuffer(lua_State *L)
{
	int channel = lua_tonumber(L,1);
	bool value = lua_toboolean(L,2);
	
	chiPie.bufferMCP3008SPI[channel]=value;
	
	return 0;
}

//########################################################## Get SPI buffer value
/** \brief Reads an SPI MCP3008 chip buffer.

\param channelNumber int Channel 0-7 to be read off the MCP3008 chip.
\param bufferPos int Buffer position (0 to 999).

\ingroup LuaPie
\author Jan*/
int chiPieGetSPIBuffer(lua_State *L)
{
	int pos = lua_tonumber(L,1);
	int value = chiPie.spi->buffer[pos];

	lua_pushnumber(L,value);
	
	return 1;
}

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

\ingroup LuaPie
\author Jan*/
int chiPieInitializeSerial(lua_State *L)
{
	int baudrate = lua_tonumber(L,1);
	
	bool success = chiPie.InitializeSerial(baudrate);
	
	lua_pushboolean(L,success);
	return 1;
}

//#################################################################### Write to serial
/** \brief Writes a message to the serial port.
\param message char Message to be sent.


\ingroup LuaPie
\author Jan*/
int chiPieSerialWrite(lua_State *L)
{
	size_t length;
	const char* message = lua_tolstring(L, 1, &length);
	
	chiPie.SerialWrite(message);
	
	return 0;
}

//#################################################################### Read from serial
/** \brief Reads a message from the serial port.

\ingroup LuaPie
\author Jan*/
int chiPieSerialRead(lua_State *L)
{
	char message[255];
	message[0]='\0';
	
	bool success = chiPie.SerialRead(message);
	lua_pushstring(L,message);
	lua_pushboolean(L,success);
	return 2;
}
