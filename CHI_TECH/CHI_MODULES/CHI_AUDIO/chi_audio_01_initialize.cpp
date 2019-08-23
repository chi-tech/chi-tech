#include "chi_audio.h"
#include "../../CHI_LIB/chi_lib.h"

#define TS_IMPLEMENTATION

bool CHI_AUDIO::LoadVorbis(char* fileName)
{
	int sampleRate;

	this->vorbisLoaded = tsLoadOGG(fileName, &sampleRate);

	CHI_LIB::StringCatEx(this->name, fileName);

	this->context = tsMakeContext(0,
								  this->frequency,
								  this->latency,
								  this->bufSeconds,
								  this->numPlayingPool);

    return true;
}


bool CHI_AUDIO::CleanAudio()
{
	// Closes Buffers
	tsShutdownContext(this->context);

	// Frees Sounds
	tsFreeSound(&this->vorbisLoaded);
}

