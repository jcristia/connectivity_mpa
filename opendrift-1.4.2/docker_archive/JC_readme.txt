in version 1.4.2, this is were the dockerfile was

I am still using 1.4.2, but I am using the dockerfile from the newest version 1.5.4.

The reason for this is that when I create my singularity image, in the old version the opendrift conda environment is a separate environment form the base env, whereas in the new dockerfile, the opendrift env is set at the base environment. When working with the singularity image, I could not figure out how to activate a different environment.