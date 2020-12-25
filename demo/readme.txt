%%%%%%%%%%%Room Information%%%%%%%%%%%%%%%%%%%%%%%%
room size: 				5 x 6 x 3 m
RT60: 					700 ms
distance between source and microphones: 	3 m
microphone spacing: 			10 cm
height of speakers:				1.6m
height of omni-directional microphone:		1 m

%%%%%%%%%%%Algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADA-RLS, RLS-based MCLP dereverberation algorithm [1];
cADA-z. constrained RLS-based MCLP dereverberation algorithm [1];
ADA-QRRLS, QRRLS-based MCLP dereverberation algorithm [2];
ADA-VFFQRRRLS, proposed VFFQRRLS-based MCLP dereverberation algorithm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sample1:   single speaker in the direction of 90°;

Sample2:   two speakers (1st speaker: 60°，0~10s;  2nd speaker: 90°，10~40s);

Note: Subfolder:
  reverberant input: unprocessed recordings (simulated by RIR generator [3]);
  reference: recordings of clean speech (reference signal, without reverberation); 
  processed results: processed speech signals by 4 adaptive MCLP dereverberation algorithms with different values of forgetting factor \gamma.

Reference:
[1] Jukić, Ante, Toon van Waterschoot, and Simon Doclo. "Adaptive speech dereverberation using constrained sparse multichannel linear prediction." IEEE Signal Processing Letters 24.1 (2016): 101-105.
[2] Sun, Xuguang, Yi Zhou, and Xiaofeng Shu. "Multi-Channel Linear Prediction Speech Dereverberation Algorithm Based on QR-RLS Adaptive Filter." Proceedings of the 3rd International Conference on Multimedia Systems and Signal Processing. 2018.
[3] Habets, Emanuel AP. "Room impulse response generator." Technische Universiteit Eindhoven, Tech. Rep 2.2.4 (2006): 1.
