#include "chi_audio.h"

//############################################################################# Default constructor
/** Default constructor.*/
CHI_AUDIO::CHI_AUDIO()
{
	this->frequency			= 44000;
	this->latency			= 15;
	this->bufSeconds		= 5;
	this->playingPool		= 1;
	this->numPlayingPool	= this->playingPool ? 5:0;

	this->context			= NULL;

	this->name[0] 			= '\0';
}


//############################################################################# Default destructor
/** Default destructor .*/
CHI_AUDIO::~CHI_AUDIO()
{

}
