#!usr/bin/python
#coding=utf-8

import re
import os
import numpy as np
import numpy.ma as ma
#from mlab.releases import R2013a as mlab
#from matlab import matlabroot

def array_to_binary_file(raw_data, output_file_name):
	tmp_data = raw_data.astype('float32')
	data = np.array(tmp_data, 'float32')
               
	fid = open(output_file_name, 'wb')
	data.tofile(fid)
	fid.close()

def extractSamples_fromUntied(neutral_untied_file,phoneme):
	content = open(neutral_untied_file,'r').readlines()
	output = open(phoneme+'_in_neutral','w')
	for i in range(len(content)):
		if '~h' in content[i]:
			full_label = (content[i].strip()).split(' ')[1][1:-1]
			A_ind = full_label.find('A:')
			line = full_label[0:A_ind]
			s_ind = line.find('-')
			e_ind = line.find('+')
			current_phone = line[s_ind+1:e_ind]
			#print current_phone
			if current_phone == phoneme:
				output.write(str(i+1)+'	'+full_label+'\n')
	output.close()

def extractHMMIdx(neutral_untied_file,hmm_index_file):
	content = open(neutral_untied_file,'r').readlines()
	output = open(hmm_index_file,'w')
	for i in range(len(content)):
		if '~h' in content[i]:
			output.write(str(i+1)+'\n')
	output.close()


def extractMatrix_fromUntied_simplified(workdir,neutral_untied_file,neutral_hmms_index_file,phoneme,mgc_dim,lf0_dim,bap_dim):

	content = open(neutral_untied_file,'r').readlines()
	indexes = open(neutral_hmms_index_file,'r').readlines()
	phoneme_samples = open('./'+phoneme+'_in_neutral','r').readlines()
	num = len(phoneme_samples)
	print "phoneme " + phoneme + " has " + str(num) + " samples"
	dim = mgc_dim*3 + lf0_dim*3 + bap_dim*3
	if os.path.exists(os.path.join(workdir,phoneme)):
		pass
	else:
		os.mkdir(os.path.join(workdir,phoneme))

	#for each state [2,3,4,5,6]
	for  i in range(2,7):
		output_state_mean = open(os.path.join(workdir,phoneme+'/neutral_'+phoneme+'_state'+str(i)+'_mean'),'w')
		output_state_variance = open(os.path.join(workdir,phoneme+'/neutral_'+phoneme+'_state'+str(i)+'_variance'),'w')

		phoneme_appearance_time = 0
		mean_state = np.zeros((num,dim))
		variance_state = np.zeros((num,dim))

		for n in range(len(phoneme_samples)):
			#print j
			phoneme_index = phoneme_samples[n].strip().split('	')[0]
			phoneme_full_label = phoneme_samples[n].strip().split('	')[1]

			for m in range(len(indexes)):
				if indexes[m].strip() == phoneme_index:
					#print phoneme_index
					hmm_start_index = int(indexes[m].strip())-1
			
					if m == len(indexes)-1:
						hmm_end_index = -1
					else:
						hmm_end_index = int(indexes[m+1].strip())-2
			
					hmm_tmp = content[hmm_start_index:hmm_end_index]

					full_label = (content[hmm_start_index].strip()).split(' ')[1][1:-1]
					print phoneme + " sample " + str(n)
					print full_label
			
					assert full_label == phoneme_full_label
					
					#for each stream [1,2,3,4,5]
					for j in range(1,6):
						#print i
						#print j
						stream_indexes = []
						mean_tmp = []
						variance_tmp = []
						for k in range(len(hmm_tmp)):
							#print k
							if '<STREAM> '+str(j) in hmm_tmp[k]:
								stream_indexes.append(k)
						#print stream_indexes
						# stream_indexes was none because 
						if j == 1 or j == 5:
							mean_tmp = hmm_tmp[stream_indexes[i-2]+2].strip().split(' ')
							variance_tmp = hmm_tmp[stream_indexes[i-2]+4].strip().split(' ')
							for t in range(len(mean_tmp)):
								mean_tmp[t] = float(mean_tmp[t])
								variance_tmp[t] = float(variance_tmp[t])
						elif j == 2 or j == 3 or j == 4:
							mean_tmp = float(hmm_tmp[stream_indexes[i-2]+4].strip())
							variance_tmp = float(hmm_tmp[stream_indexes[i-2]+6].strip())

						if j == 1:
							mean_state[phoneme_appearance_time,0:mgc_dim*3] = mean_tmp
							variance_state[phoneme_appearance_time,0:mgc_dim*3] = variance_tmp
						elif j == 5:
							mean_state[phoneme_appearance_time,mgc_dim*3+lf0_dim*3:dim] = mean_tmp
							variance_state[phoneme_appearance_time,mgc_dim*3+lf0_dim*3:dim] = variance_tmp
						else:
							mean_state[phoneme_appearance_time,mgc_dim*3+lf0_dim*(j-2):mgc_dim*3+lf0_dim*(j-2)+1] = mean_tmp
							variance_state[phoneme_appearance_time,mgc_dim*3+lf0_dim*(j-2):mgc_dim*3+lf0_dim*(j-2)+1] = variance_tmp

			#if phoneme_appearance_time > 1:
			#	assert variance_state[phoneme_appearance_time,:] != variance_state[phoneme_appearance_time-1]
			
			phoneme_appearance_time = phoneme_appearance_time+1
			print "phoneme_appearance_time "+ str(phoneme_appearance_time)

		np.savetxt(output_state_mean, mean_state)
		np.savetxt(output_state_variance, variance_state)

		output_state_mean.close()
		output_state_variance.close()



def extractMatrix_fromReclustered_simplified(workdir,happy_reclustered_file,happy_hmms_index_file,phoneme,mgc_dim,lf0_dim,bap_dim):
	lf0_map = {}
	mgc_map = {}
	bap_map = {}

	content = open(happy_reclustered_file,'r').readlines()
	#lf0_map_output = open('lf0_map','w')
	#mgc_map_output = open('mgc_map','w')
	#bap_map_output = open('bap_map','w')
	for i in range(len(content)):
		if '~p' in content[i] and 'lf0' in content[i].strip().split(' ')[1][1:-1] and '<MIXTURE>' in content[i+3]:
			key = content[i].strip().split(' ')[1][1:-1]
			#print key
			lf0_map[key] = i+1
			#print i+1
			#lf0_map_output.write(key+'	'+str(lf0_map[key])+'\n')
		if '~p' in content[i] and 'mgc' in content[i].strip().split(' ')[1][1:-1] and 'GCONST' in content[i-1]:
			key = content[i].strip().split(' ')[1][1:-1]
			#print key
			mgc_map[key] = i+1
			#print i+1
			#mgc_map_output.write(key+'	'+str(mgc_map[key])+'\n')
		if '~p' in content[i] and 'bap' in content[i].strip().split(' ')[1][1:-1] and 'GCONST' in content[i-1]:
			key = content[i].strip().split(' ')[1][1:-1]
			#print key
			bap_map[key] = i+1
			#print i+1
			#bap_map_output.write(key+'	'+str(bap_map[key])+'\n')

	phoneme_samples = open(phoneme+'_in_neutral','r').readlines()
	indexes = open(happy_hmms_index_file,'r').readlines()
	num = len(phoneme_samples)
	#print num
	dim = mgc_dim*3 + lf0_dim*3 + bap_dim*3

	#for each state [2,3,4,5,6]
	for i in range(2,7):
		output_state_mean = open(os.path.join(workdir,phoneme+'/happy_'+phoneme+'_state'+str(i)+'_mean'),'w')
		output_state_variance = open(os.path.join(workdir,phoneme+'/happy_'+phoneme+'_state'+str(i)+'_variance'),'w')

		phoneme_appearance_time = 0

		mean_state = np.zeros((num,dim))
		variance_state = np.zeros((num,dim))
	
		for n in range(len(phoneme_samples)):
			#phoneme_index = phoneme_samples[n].strip().split('	')[0]
			phoneme_full_label = phoneme_samples[n].strip().split('	')[1]
			for m in range(len(indexes)):
				#if indexes[m] == phoneme_index:
				full_label = content[int(indexes[m])-1].strip().split(' ')[1][1:-1]
				#print full_label
				if full_label == phoneme_full_label:

					#print "we found the sample"+str(phoneme_appearance_time+1)
					#print i
					
					state_index = int(indexes[m])-1 + 3 + 12*(i-2)
					#print int(indexes[m])-1
					#print state_index
					#print content[state_index]


					for j in range(1,6):
						state_stream = content[state_index+1+2*j].strip().split(' ')[1][1:-1]
						#print state_stream
						if j == 1:
							mgc_index = mgc_map[state_stream]
							mgc_mean = content[mgc_index-1+3].strip().split(' ')
							mgc_variance = content[mgc_index-1+5].strip().split(' ')
							for k in range(len(mgc_mean)):
								mgc_mean[k] = float(mgc_mean[k])
								mgc_variance[k] = float(mgc_variance[k])
							mean_state[phoneme_appearance_time,0:mgc_dim*3] = mgc_mean
							variance_state[phoneme_appearance_time,0:mgc_dim*3] = mgc_variance
						elif j == 2 or j == 3 or j == 4:
							lf0_index = lf0_map[state_stream]
							lf0_mean = float(content[lf0_index-1+5].strip())
							lf0_variance = float(content[lf0_index-1+7].strip())
							mean_state[phoneme_appearance_time,mgc_dim*3+lf0_dim*(j-2):mgc_dim*3+lf0_dim*(j-2)+1] = lf0_mean
							variance_state[phoneme_appearance_time,mgc_dim*3+lf0_dim*(j-2):mgc_dim*3+lf0_dim*(j-2)+1] = lf0_variance
						else:
							bap_index = bap_map[state_stream]
							bap_mean = content[bap_index-1+3].strip().split(' ')
							bap_variance = content[bap_index-1+5].strip().split(' ')
							for t in range(len(bap_mean)):
								bap_mean[t] = float(bap_mean[t])
								bap_variance[t] = float(bap_variance[t])
							mean_state[phoneme_appearance_time,mgc_dim*3+lf0_dim*3:dim] = bap_mean
							variance_state[phoneme_appearance_time,mgc_dim*3+lf0_dim*3:dim] = bap_variance							




			phoneme_appearance_time += 1

		np.savetxt(output_state_mean, mean_state)
		np.savetxt(output_state_variance, variance_state)

		output_state_mean.close()
		output_state_variance.close()

def train_convert_matrix(dir,phoneme,mgc_dim,lf0_dim,bap_dim):
	dim = mgc_dim*3 + lf0_dim*3 + bap_dim*3;
	for i in range(2,7):
		mean_tran = np.zeros((dim,dim))
		variance_tran = np.zeros((dim,dim))

		neutral_mean_matrix = np.loadtxt(os.path.join(dir,'neutral_'+phoneme+'_state'+str(i)+'_mean'))
		happy_mean_matrix = np.loadtxt(os.path.join(dir,'happy_'+phoneme+'_state'+str(i)+'_mean'))
		x_mean = np.array(neutral_mean_matrix, 'float32')
		y_mean = np.array(happy_mean_matrix, 'float32')
		print x_mean.shape
		print y_mean.shape
	
		for j in range(dim):
			mean_tran[:,j]=np.dot(np.dot(np.linalg.pinv(np.dot(x_mean.T,x_mean)),x_mean.T),y_mean[:,j]);

		fid_mean = open(os.path.join(dir,'neutral2happy_'+phoneme+'_state'+str(i)+'_mean_transformation'),'w');
		#lines
		n = mean_tran.shape[0]
		#columns
		m = mean_tran.shape[1]
		if n != m:
			print('mean_tran is not square matrix')
		
		np.savetxt(fid_mean, mean_tran)
		fid_mean.close()

		neutral_variance_matrix = np.loadtxt(os.path.join(dir,'neutral_'+phoneme+'_state'+str(i)+'_variance'))
		happy_variance_matrix = np.loadtxt(os.path.join(dir,'happy_'+phoneme+'_state'+str(i)+'_variance'))

		x_variance = np.array(neutral_variance_matrix, 'float32')
		y_variance = np.array(happy_variance_matrix, 'float32')

		#mask_x=x_variance_ori<1e-6
		#mask_y=y_variance_ori<1e-6
		#masked_x_variance = ma.array(x_variance_ori,mask=mask_x,fill_value=1e-4)
		#masked_y_variance = ma.array(y_variance_ori,mask=mask_y,fill_value=1e-4)
		#x_variance = masked_x_variance.filled()*10000
		#y_variance = masked_y_variance.filled()*10000
		#print min(x_variance)
		#print max(x_variance)

		#exit()

		print x_variance.shape
		print y_variance.shape
	
		for j in range(dim):
			variance_tran[:,j]=np.dot(np.dot(np.linalg.pinv(np.dot(x_variance.T,x_variance)),x_variance.T),y_variance[:,j]);

		fid_variance = open(os.path.join(dir,'neutral2happy_'+phoneme+'_state'+str(i)+'_variance_transformation'),'w');
		#lines
		n = variance_tran.shape[0]
		#columns
		m = variance_tran.shape[1]
		if n != m:
			print('variance_tran is not square matrix')
		
		np.savetxt(fid_variance, variance_tran)
		fid_variance.close()


def convert(workdir,full_label_name,neutral_reclustered_file,neutral_untied_file,happy_reclustered_file,neutral_hmms_index_file,happy_hmms_index_file,mgc_dim,lf0_dim,bap_dim):
	lf0_map = {}
	mgc_map = {}
	bap_map = {}
	# use the result of re_clustered_all.mmf.1mix, while we use the results of untied.mmf to train the convert matrix
	content = open(neutral_reclustered_file,'r').readlines()
	full_label = open(full_label_name,'r').readlines()
	file_name = full_label_name.split('/')[-1][0:-4]
	print "we need to convert this file: "+file_name
	output_ori_neutral_variance = open(os.path.join(workdir,file_name+'_ori_neutral_state_variance'),'w')
	output_ori_neutral_mean = open(os.path.join(workdir,file_name+'_ori_neutral_state_mean'),'w')
	output_pre_happy_variance = open(os.path.join(workdir,file_name+'_pre_happy_state_variance'),'w')
	output_pre_happy_mean = open(os.path.join(workdir,file_name+'_pre_happy_state_mean'),'w')

	for i in range(len(content)):
		if '~p' in content[i] and 'lf0' in content[i].strip().split(' ')[1][1:-1] and 'GCONST' in content[i-1]:
			key = content[i].strip().split(' ')[1][1:-1]
			lf0_map[key] = i+1

		if '~p' in content[i] and 'mgc' in content[i].strip().split(' ')[1][1:-1] and 'GCONST' in content[i-1]:
			key = content[i].strip().split(' ')[1][1:-1]
			mgc_map[key] = i+1

		if '~p' in content[i] and 'bap' in content[i].strip().split(' ')[1][1:-1] and 'GCONST' in content[i-1]:
			key = content[i].strip().split(' ')[1][1:-1]
			bap_map[key] = i+1

	
	dim = mgc_dim*3 + lf0_dim*3 + bap_dim*3
	mean = np.zeros((len(full_label)*5,dim))
	variance = np.zeros((len(full_label)*5,dim))
	pre_mean = np.zeros((len(full_label)*5,dim))
	pre_variance = np.zeros((len(full_label)*5,dim))
	index = 0
	#for each line in full label, we need to collect 10 matrixes in both neutral and happy mmf(mean and variance for 5 states), and train 10 transform matrixes
	for line in full_label:
		A_ind = (line.strip()).find('A:')
		line_part = (line.strip())[0:A_ind]
		s_ind = line_part.find('-')
		e_ind = line_part.find('+')
		phoneme = line_part[s_ind+1:e_ind]
		print "the current phoneme is: "+phoneme
		#print os.path.join(workdir,phoneme+'in_neutral')
		if os.path.exists(os.path.join(workdir,phoneme+'_in_neutral')) and phoneme != 'pau':
			print phoneme+'_in_neutral is there!'
		elif phoneme == 'pau':
			print 'we do not convert pau'
		else:
			print 'we need to make it'
			#collect all the samples of the specific phoneme in the neutral untied.mmf
			extractSamples_fromUntied(neutral_untied_file,phoneme)
			#collect the 10 matrixes (mean and variance matrixes of 5 states) for the specific phoneme in the neutral untied.mmf 
			extractMatrix_fromUntied_simplified(workdir,neutral_untied_file,neutral_hmms_index_file, phoneme,mgc_dim,lf0_dim,bap_dim)
			#collect the 10 matrixes (mean and variance matrixes of 5 states) for the specific phoneme in the happy re_clustered_all.mmf.1mix
			#print "we've been here"
			extractMatrix_fromReclustered_simplified(workdir,happy_reclustered_file,happy_hmms_index_file,phoneme,mgc_dim,lf0_dim,bap_dim)
			#now, there is a folder named by the specific phoneme, and in the folder, there are 20 matrixes waiting to train the convert transform matrixes.
			train_convert_matrix(os.path.join(workdir,phoneme),phoneme,mgc_dim,lf0_dim,bap_dim)

		#each phoneme has its own transformation matrix, and each state in one phoneme has its own transformation too.
		#but the convert matrixes of the five states belonging to one phoneme are in the same folder, so we use set transformation of one phoneme as a block
		mean_tmp = np.zeros((5,dim))
		variance_tmp = np.zeros((5,dim))
		pre_mean_tmp = np.zeros((5,dim))
		pre_variance_tmp = np.zeros((5,dim))
		#print np.shape(mean_tmp)
		for i in range(len(content)):
			#print i
			if '~h' in content[i]:
				if (content[i].strip()).split(' ')[1][1:-1] == line.strip():
					print "the corresponding full_label is in line: "+str(i)
					print phoneme+" can be found in neutral reclustered file successfully"
					#for each state [2,3,4,5,6]
					for n in range(2,7):
						state_index = i + 3 + 12*(n-2)
						#print content[state_index]
						#for each stream [1,2,3,4,5]
						for m in range(1,6):
							state_stream = content[state_index+1+2*m].strip().split(' ')[1][1:-1]
							#print state_stream
							if m == 1:
								mgc_index = mgc_map[state_stream]
								mgc_mean = content[mgc_index-1+3].strip().split(' ')
								mgc_variance = content[mgc_index-1+5].strip().split(' ')
								for k in range(len(mgc_mean)):
									mgc_mean[k] = float(mgc_mean[k])
									mgc_variance[k] = float(mgc_variance[k])
								mean_tmp[n-2,0:mgc_dim*3] = mgc_mean
								variance_tmp[n-2,0:mgc_dim*3] = mgc_variance
							elif m == 2 or m == 3 or m == 4:
								lf0_index = lf0_map[state_stream]
								#here we just use the first mixture of lf0 model (even though the weight of the first mixture is not 1.000000e+1)
								lf0_mean = float(content[lf0_index-1+5].strip())
								lf0_variance = float(content[lf0_index-1+7].strip())
								mean_tmp[n-2,mgc_dim*3+lf0_dim*(m-2):mgc_dim*3+lf0_dim*(m-2)+1] = lf0_mean
								variance_tmp[n-2,mgc_dim*3+lf0_dim*(m-2):mgc_dim*3+lf0_dim*(m-2)+1] = lf0_variance
							else:
								bap_index = bap_map[state_stream]
								bap_mean = content[bap_index-1+3].strip().split(' ')
								bap_variance = content[bap_index-1+5].strip().split(' ')
								for k in range(len(bap_mean)):
									bap_mean[k] = float(bap_mean[k])
									bap_variance[k] = float(bap_variance[k])
								mean_tmp[n-2,mgc_dim*3+lf0_dim*3:dim] = bap_mean
								variance_tmp[n-2,mgc_dim*3+lf0_dim*3:dim] = bap_variance
						#if current state belongs to pau, then we do nothing
						#if current state doesn't belongs to pau, then we do conversion
						#variance_tmp[n-2,:] = variance_tmp[n-2,:]*10000
						if phoneme == 'pau':
							pass
						else:
							mean_tran_matrix = np.loadtxt(os.path.join(workdir,phoneme+'/neutral2happy_'+phoneme+'_state'+str(n)+'_mean_transformation'))
							variance_tran_matrix = np.loadtxt(os.path.join(workdir,phoneme+'/neutral2happy_'+phoneme+'_state'+str(n)+'_mean_transformation'))
							#print np.shape(mean_tran_matrix)
							#print mean_tmp[n-2,:]
							#print np.dot(np.ones((1,93)),np.ones((93,1)))
							for d in range(0,dim):
								pre_mean_tmp[n-2,d] = np.dot(mean_tmp[n-2,:],mean_tran_matrix[:,d])
								pre_variance_tmp[n-2,d] = np.dot(variance_tmp[n-2,:],variance_tran_matrix[:,d])

					#print np.shape(mean_tmp)
					#print np.shape(variance_tmp)
					if phoneme == 'pau':
						mean[index*5:(index+1)*5,:] = mean_tmp
						variance[index*5:(index+1)*5,:] = variance_tmp
						pre_mean[index*5:(index+1)*5,:] = mean_tmp
						pre_variance[index*5:(index+1)*5,:] = variance_tmp
					else:
						mean[index*5:(index+1)*5,:] = mean_tmp
						variance[index*5:(index+1)*5,:] = variance_tmp
						pre_mean[index*5:(index+1)*5,:] = pre_mean_tmp
						pre_variance[index*5:(index+1)*5,:] = pre_variance_tmp

		index += 1
		print "the "+str(index)+" line in full label has been converted"

	#print np.shape(mean)
	#print np.shape(variance)
	np.savetxt(output_ori_neutral_mean, mean)
	np.savetxt(output_ori_neutral_variance, variance)
	np.savetxt(output_pre_happy_mean, pre_mean)
	np.savetxt(output_pre_happy_variance, pre_variance)

	output_ori_neutral_mean.close()
	output_ori_neutral_variance.close()
	output_pre_happy_mean.close()
	output_pre_happy_variance.close()


#I've changed it for a lot of times when debugging
def addDurInfo(workdir, full_label_name, dur_reclustered_file, pre_happy_state_mean, pre_happy_state_variance,mgc_dim,lf0_dim,bap_dim):
	
	file_name = full_label_name.split('/')[-1][0:-4]
	tmp = (pre_happy_state_mean.split('/')[-1]).split('_')[1]+'_'+(pre_happy_state_mean.split('/')[-1]).split('_')[2]
	dim = mgc_dim*3+lf0_dim*3+bap_dim*3

	
	output_pre_happy_mgc_mean_variance_frame_level_name = os.path.join(workdir,file_name+'_'+tmp+'_frame_mean_variance_mgc')
	output_pre_happy_lf0_mean_variance_frame_level_name = os.path.join(workdir,file_name+'_'+tmp+'_frame_mean_variance_lf0')
	output_pre_happy_bap_mean_variance_frame_level_name = os.path.join(workdir,file_name+'_'+tmp+'_frame_mean_variance_bap')

	dur_map={}
	content = open(dur_reclustered_file,'r').readlines()
	for i in range(len(content)):
		if '~s' in content[i] and '<STREAM>' in content[i+1]:
			key = content[i].strip().split(' ')[1][1:-1]
			#when we collect lf0_map or bap_map or mgc_bap, the index always add 1, that's because i wanna check in the document, ahha
			dur_map[key] = i

	full_label = open(full_label_name,'r').readlines()
	mean = np.loadtxt(pre_happy_state_mean)
	variance = np.loadtxt(pre_happy_state_variance)

	index = 0
	repeat_num_total = 0
	state_dur_map={}
	for line in full_label:
		for i in range(len(content)):
			if '~h' in content[i] and content[i].strip().split(' ')[1][1:-1] == line.strip():
				print "I found it! "+ str(index+1)
				dur_info = content[i+4].strip().split(' ')[1][1:-1]
				dur_index = dur_map[dur_info]
				#print dur_index
				for j in range(2,7):
					dur_mean = content[int(dur_index)+3+6*(j-2)]
					dur_variance = content[int(dur_index)+5+6*(j-2)]
					repeat_num = int(round(float(dur_mean)))
					#print "state "+str(j)
					#print dur_mean
					#print repeat_num
					repeat_num_total = repeat_num_total + repeat_num
					key = index*5+(j-2)
					#print key
					state_dur_map[key] = repeat_num

		index += 1
	print repeat_num_total

	final_mean = np.zeros((repeat_num_total,mean.shape[1]))
	final_variance = np.zeros((repeat_num_total,variance.shape[1]))
	final_mgc_mean_variance = np.zeros((repeat_num_total,mgc_dim*3*2))
	final_lf0_mean_variance = np.zeros((repeat_num_total,lf0_dim*3*2))
	final_bap_mean_variance = np.zeros((repeat_num_total,bap_dim*3*2))
	current_repeat_num = 0

	for n in range(mean.shape[0]):
		for k in range(state_dur_map[n]):
			final_mean[current_repeat_num+k,:]=mean[n,:]
			final_variance[current_repeat_num+k,:]=variance[n,:]
		current_repeat_num = current_repeat_num + state_dur_map[n]

	final_mgc_mean_variance[:,0:mgc_dim*3]=final_mean[:,0:mgc_dim*3]
	final_mgc_mean_variance[:,mgc_dim*3:mgc_dim*3*2]=final_variance[:,0:mgc_dim*3]
	#final_mgc_mean_variance[:,mgc_dim*3:mgc_dim*3+mgc_dim]=np.abs(final_variance[:,0:mgc_dim])
	#final_mgc_mean_variance[:,mgc_dim*3+mgc_dim:mgc_dim*3*2]=final_variance[:,mgc_dim:mgc_dim*3]
	final_lf0_mean_variance[:,0:lf0_dim*3]=final_mean[:,mgc_dim*3:mgc_dim*3+lf0_dim*3]
	final_lf0_mean_variance[:,lf0_dim*3:lf0_dim*3*2]=final_variance[:,mgc_dim*3:mgc_dim*3+lf0_dim*3]
	final_bap_mean_variance[:,0:bap_dim*3]=final_mean[:,mgc_dim*3+lf0_dim*3:dim]
	final_bap_mean_variance[:,bap_dim*3:bap_dim*3*2]=final_variance[:,mgc_dim*3+lf0_dim*3:dim]

	array_to_binary_file(final_mgc_mean_variance, output_pre_happy_mgc_mean_variance_frame_level_name)
	array_to_binary_file(final_lf0_mean_variance, output_pre_happy_lf0_mean_variance_frame_level_name)
	array_to_binary_file(final_bap_mean_variance, output_pre_happy_bap_mean_variance_frame_level_name)

	output_pre_happy_mgc_mean_variance_frame_level = open(output_pre_happy_mgc_mean_variance_frame_level_name+'.txt','w')
	output_pre_happy_lf0_mean_variance_frame_level = open(output_pre_happy_lf0_mean_variance_frame_level_name+'.txt','w')
	output_pre_happy_bap_mean_variance_frame_level = open(output_pre_happy_bap_mean_variance_frame_level_name+'.txt','w')
	output_pre_happy_mean_frame_level = open(os.path.join(workdir,file_name+'_ori_neutral_frame_mean'),'w')
	output_pre_happy_variance_frame_level = open(os.path.join(workdir,file_name+'_ori_neutral_frame_variance'),'w')

	np.savetxt(output_pre_happy_mgc_mean_variance_frame_level, final_mgc_mean_variance)
	np.savetxt(output_pre_happy_lf0_mean_variance_frame_level, final_lf0_mean_variance)
	np.savetxt(output_pre_happy_bap_mean_variance_frame_level, final_bap_mean_variance)
	np.savetxt(output_pre_happy_mean_frame_level, final_mean)
	np.savetxt(output_pre_happy_variance_frame_level, final_variance)

	output_pre_happy_mgc_mean_variance_frame_level.close()
	output_pre_happy_lf0_mean_variance_frame_level.close()
	output_pre_happy_bap_mean_variance_frame_level.close()
	output_pre_happy_mean_frame_level.close()
	output_pre_happy_variance_frame_level.close()

def convert_check(workdir,full_label_name,neutral_reclustered_file,neutral_untied_file,happy_reclustered_file,mgc_dim,lf0_dim,bap_dim):
	lf0_map = {}
	mgc_map = {}
	bap_map = {}
	# use the result of re_clustered_all.mmf.1mix, while we use the results of untied.mmf to train the convert matrix
	content = open(happy_reclustered_file,'r').readlines()
	full_label = open(full_label_name,'r').readlines()
	file_name = full_label_name.split('/')[-1][0:-4]
	print "we need to collect this file: "+file_name
	output_ori_happy_variance = open(os.path.join(workdir,file_name+'_ori_happy_state_variance'),'w')
	output_ori_happy_mean = open(os.path.join(workdir,file_name+'_ori_happy_state_mean'),'w')
	#output_pre_happy_variance = open(os.path.join(workdir,file_name+'_pre_happy_state_variance'),'w')
	#output_pre_happy_mean = open(os.path.join(workdir,file_name+'_pre_happy_state_mean'),'w')

	for i in range(len(content)):
		if '~p' in content[i] and 'lf0' in content[i].strip().split(' ')[1][1:-1] and 'GCONST' in content[i-1]:
			key = content[i].strip().split(' ')[1][1:-1]
			lf0_map[key] = i+1

		if '~p' in content[i] and 'mgc' in content[i].strip().split(' ')[1][1:-1] and 'GCONST' in content[i-1]:
			key = content[i].strip().split(' ')[1][1:-1]
			mgc_map[key] = i+1

		if '~p' in content[i] and 'bap' in content[i].strip().split(' ')[1][1:-1] and 'GCONST' in content[i-1]:
			key = content[i].strip().split(' ')[1][1:-1]
			bap_map[key] = i+1

	
	dim = mgc_dim*3 + lf0_dim*3 + bap_dim*3
	mean = np.zeros((len(full_label)*5,dim))
	variance = np.zeros((len(full_label)*5,dim))
	#pre_mean = np.zeros((len(full_label)*5,dim))
	#pre_variance = np.zeros((len(full_label)*5,dim))
	index = 0
	#for each line in full label, we need to collect 10 matrixes in both neutral and happy mmf(mean and variance for 5 states), and train 10 transform matrixes
	for line in full_label:
		A_ind = (line.strip()).find('A:')
		line_part = (line.strip())[0:A_ind]
		s_ind = line_part.find('-')
		e_ind = line_part.find('+')
		phoneme = line_part[s_ind+1:e_ind]
		print "the current phoneme is: "+phoneme
		#print os.path.join(workdir,phoneme+'in_neutral')

		#each phoneme has its own transformation matrix, and each state in one phoneme has its own transformation too.
		#but the convert matrixes of the five states belonging to one phoneme are in the same folder, so we use set transformation of one phoneme as a block
		mean_tmp = np.zeros((5,dim))
		variance_tmp = np.zeros((5,dim))
		#pre_mean_tmp = np.zeros((5,dim))
		#pre_variance_tmp = np.zeros((5,dim))
		#print np.shape(mean_tmp)
		for i in range(len(content)):
			#print i
			if '~h' in content[i]:
				if (content[i].strip()).split(' ')[1][1:-1] == line.strip():
					print "the corresponding full_label is in line: "+str(i)
					print phoneme+" can be found in happy reclustered file successfully"
					#for each state [2,3,4,5,6]
					for n in range(2,7):
						state_index = i + 3 + 12*(n-2)
						#print content[state_index]
						#for each stream [1,2,3,4,5]
						for m in range(1,6):
							state_stream = content[state_index+1+2*m].strip().split(' ')[1][1:-1]
							#print state_stream
							if m == 1:
								mgc_index = mgc_map[state_stream]
								mgc_mean = content[mgc_index-1+3].strip().split(' ')
								mgc_variance = content[mgc_index-1+5].strip().split(' ')
								for k in range(len(mgc_mean)):
									mgc_mean[k] = float(mgc_mean[k])
									mgc_variance[k] = float(mgc_variance[k])
								mean_tmp[n-2,0:mgc_dim*3] = mgc_mean
								variance_tmp[n-2,0:mgc_dim*3] = mgc_variance
							elif m == 2 or m == 3 or m == 4:
								lf0_index = lf0_map[state_stream]
								#here we just use the first mixture of lf0 model (even though the weight of the first mixture is not 1.000000e+1)
								lf0_mean = float(content[lf0_index-1+5].strip())
								lf0_variance = float(content[lf0_index-1+7].strip())
								mean_tmp[n-2,mgc_dim*3+lf0_dim*(m-2):mgc_dim*3+lf0_dim*(m-2)+1] = lf0_mean
								variance_tmp[n-2,mgc_dim*3+lf0_dim*(m-2):mgc_dim*3+lf0_dim*(m-2)+1] = lf0_variance
							else:
								bap_index = bap_map[state_stream]
								bap_mean = content[bap_index-1+3].strip().split(' ')
								bap_variance = content[bap_index-1+5].strip().split(' ')
								for k in range(len(bap_mean)):
									bap_mean[k] = float(bap_mean[k])
									bap_variance[k] = float(bap_variance[k])
								mean_tmp[n-2,mgc_dim*3+lf0_dim*3:dim] = bap_mean
								variance_tmp[n-2,mgc_dim*3+lf0_dim*3:dim] = bap_variance

					mean[index*5:(index+1)*5,:] = mean_tmp
					variance[index*5:(index+1)*5,:] = variance_tmp


		index += 1
		print "the "+str(index)+" line in full label has been converted"

	np.savetxt(output_ori_happy_mean, mean)
	np.savetxt(output_ori_happy_variance, variance)

	output_ori_happy_mean.close()
	output_ori_happy_variance.close()


if __name__ == '__main__':
	#check if untid.mmf and fullcontext.mmf have the same number 
	#extractSamples_fromUntied('/home/will/workspace/hts/HTS-neutral_wangzhe/models/qst001/ver1/cmp/fullcontext.mmf.txt', 'CNx')
	#extractHMMIdx('/home/will/workspace/hts/HTS-happy_wangzhe/models/qst001/ver1/cmp/re_clustered_all.mmf.1mix.txt','happy_reclustered_hmms_index')
	#extractHMMIdx('/home/will/workspace/hts/HTS-neutral_wangzhe/models/qst001/ver1/cmp/fullcontext.mmf.txt','neutral_fullcontext_hmms_index_new')
	#extractMatrix_fromUntied_simplified('./HTS-neutral_wangzhe/models/qst001/ver1/cmp/untied.mmf.txt', 'CNuw',25,1,5)
	  #extractMatrix_fromUntied('./HTS-neutral_wangzhe/models/qst001/ver1/cmp/untied.mmf.txt')
	  #extractMatrix_fromReclustered('./HTS-happy_wangzhe/models/qst001/ver1/cmp/re_clustered_all.mmf.1mix.txt')
	#extractMatrix_fromReclustered_simplified('./HTS-happy_wangzhe/models/qst001/ver1/cmp/re_clustered_all.mmf.1mix.txt','CNuw',25,1,5)
	#train_convert_matrix('/home/will/workspace/hts/conversion/CNx','CNx',25,1,5)
	workdir = '/home/will/workspace/hts/conversion'
	neutral_reclustered_file = '/home/will/workspace/hts/HTS-neutral_wangzhe/models/qst001/ver1/cmp/re_clustered_all.mmf.1mix.txt'
	#neutral_untied_file = '/home/will/workspace/hts/HTS-neutral_wangzhe/models/qst001/ver1/cmp/untied.mmf.txt'
	neutral_untied_file = '/home/will/workspace/hts/HTS-neutral_wangzhe/models/qst001/ver1/cmp/fullcontext.mmf.txt'
	neutral_hmms_index_file = '/home/will/workspace/hts/conversion/neutral_fullcontext_hmms_index'
	happy_hmms_index_file = '/home/will/workspace/hts/conversion/happy_reclustered_hmms_index'
	happy_reclustered_file = '/home/will/workspace/hts/HTS-happy_wangzhe/models/qst001/ver1/cmp/re_clustered_all.mmf.1mix.txt'
	full_label_name = '/home/will/workspace/hts/HTS-neutral_wangzhe/data/labels/gen/206.lab'
	dur_reclustered_file = '/home/will/workspace/hts/HTS-happy_wangzhe/models/qst001/ver1/dur/re_clustered_all.mmf.1mix.txt'
	pre_happy_state_mean = '/home/will/workspace/hts/conversion/206_pre_happy_state_mean'
	pre_happy_state_variance = '/home/will/workspace/hts/conversion/206_pre_happy_state_variance'
	ori_neutral_state_mean = '/home/will/workspace/hts/conversion/206_ori_neutral_state_mean'
	ori_neutral_state_variance = '/home/will/workspace/hts/conversion/206_ori_neutral_state_variance'
	convert(workdir,full_label_name,neutral_reclustered_file,neutral_untied_file,happy_reclustered_file,neutral_hmms_index_file,happy_hmms_index_file,25,1,5)
	#addDurInfo(workdir,full_label_name,dur_reclustered_file,ori_neutral_state_mean,ori_neutral_state_variance,25,1,5)
	#addDurInfo(workdir,full_label_name,dur_reclustered_file,pre_happy_state_mean,pre_happy_state_variance,25,1,5)
	#convert_check(workdir,full_label_name,neutral_reclustered_file,neutral_untied_file,happy_reclustered_file,25,1,5)
