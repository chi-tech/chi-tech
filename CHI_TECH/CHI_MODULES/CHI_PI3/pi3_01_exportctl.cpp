#include "chi_pi3.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


//######################################################### Export a pin
/**Exports a pin to sysfs.
\param pinNum The pin number (1-40) to be exported.
\author Jan*/
bool CHI_PI3::ExportPin(int pinNum)
{
	if ((pinNum<=0) || (pinNum>=41)){return false;}

	std::ofstream ofs;

	ofs.open("/sys/class/gpio/export");

	if (ofs.fail())
	{
		printf("Unable to export GPIO pin #%d\n",pinNum);
	}

	ofs << this->pinNumber[pinNum];

	this->pinExported[pinNum]	= true;

	ofs.close();

	return true;
}

//######################################################### Export a pin
/**UnExports a pin to sysfs.
\param pinNum The pin number (1-40) to be exported.
\author Jan*/
bool CHI_PI3::UnExportPin(int pinNum)
{
	if ((pinNum<=0) || (pinNum>=41)){return false;}

	std::ofstream ofs;

	ofs.open("/sys/class/gpio/unexport");

	if (ofs.fail())
	{
		printf("Unable to export GPIO pin #%d\n",pinNum);
	}

	ofs << this->pinNumber[pinNum];

	this->pinExported[pinNum]	= false;

	ofs.close();

	return true;
}

//######################################################### Set pin mode
/**Sets a pin's mode to either read or write.
\param pinNum The pin number (1-40) to be exported.
\param mode    PI3_READ for an input pin or, PI3_WRITE for an output pin.

\note
Writing to pins designated as READ has no effect. Reading from pins designated as
WRITE will give the current value of the pin (i.e. WRITE does not work on READ-pins but READ
works on WRITE-pins).

\author Jan*/
bool CHI_PI3::SetPinMode(int pinNum, int mode)
{
	if (!this->pinExported[pinNum]){this->ExportPin(pinNum);}

	if ((pinNum<=0) || (pinNum>=41)){return false;}

	std::string setdir_str ="/sys/class/gpio/gpio" + this->pinNumber[pinNum] + "/direction";

	std::ofstream ofs;

	ofs.open(setdir_str.c_str());

	if (ofs.fail())
	{
		printf("Unable to set direction of GPIO pin #%d\n",pinNum);
	}
	if (mode==PI3_READ)
	{
		ofs << "in";
	}
	else if (mode==PI3_WRITE)
	{
		ofs << "out";
	}
	
	ofs.close();
	return true;
}

//######################################################### Set pin value
/**Sets a pin's value to either high or low.
\param pinNum The pin number (1-40) to be exported.
\param value     PI3_LOW for 0 or, PI3_HIGH for 1.
\author Jan*/
bool CHI_PI3::SetPinValue(int pinNum, int value)
{
	if ((pinNum<=0) || (pinNum>=41)){return false;}

	std::string setval_str ="/sys/class/gpio/gpio" + this->pinNumber[pinNum] + "/value";

	std::ofstream ofs;

	ofs.open(setval_str.c_str());

	if (ofs.fail())
	{
		printf("Unable to set value of GPIO pin #%d\n",pinNum);
	}
	if (value==PI3_LOW)
	{
		ofs << "0";
	}
	else if (value==PI3_HIGH)
	{
		ofs << "1";
	}
	
	ofs.close();
	return true;
}

//######################################################### Get pin value
/**Gets a pin's value.
\param pinNum The pin number (1-40) to be exported.
\return Either PI3_HIGH or PI3_LOW.
\author Jan*/
int CHI_PI3::GetPinValue(int pinNum)
{
	if ((pinNum<=0) || (pinNum>=41)){return false;}

	std::string val;

	std::string getval_str  ="/sys/class/gpio/gpio" + this->pinNumber[pinNum] + "/value";

	std::ifstream ofs;

	ofs.open(getval_str.c_str());
	if (ofs.fail())
	{
		printf("Unable to get value of GPIO pin #%d\n",pinNum);
	}

	ofs >> val;

	if (val=="0")
	{
		return PI3_LOW;
	}
	else 
	{
		return PI3_HIGH;
	}
	
	ofs.close();
	return 1;
}
