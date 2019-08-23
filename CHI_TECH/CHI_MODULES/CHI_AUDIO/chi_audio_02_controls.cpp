#include "chi_audio.h"



void CHI_AUDIO::PlayAudio(int loopAmount)
{
	if(&this->vorbisLoaded == NULL)
	{
		printf("Error sound file %s could not be played. \n", this->name);
	}
	else
	{
		this->vorbisSound = tsMakeDef(&this->vorbisLoaded);

		if(loopAmount > 0)
		{
			this->vorbisSound.looped = loopAmount;
			tsPlaySound(this->context, this->vorbisSound);
		}
		else
		{
			tsPlaySound(this->context, this->vorbisSound);
		}
	}
}



void CHI_AUDIO::StopAudio()
{

}



void CHI_AUDIO::VolumeAudio(float volume)
{
	this->vorbisSound.volume_left = volume;
	this->vorbisSound.volume_right = volume;
}
