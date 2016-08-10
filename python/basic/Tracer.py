#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'

_tracer_mute = 0
_tracer_error = 1
_tracer_warning = 2
_tracer_info = 3
_tracer_debug = 4
_tracer_trace = 5

level_strings={'mute':_tracer_mute,
	      'error':_tracer_error,
				'warning':_tracer_warning,
				'info':_tracer_info,
				'debug':_tracer_debug,
				'trace':_tracer_trace}

import sys

settings={}
tracers={}
initialized=False
warned_on_uninitialized=False

def init_tracer_from_cmdline(cmdline):
	global settings
	global initialized
	initialized=True

	if not cmdline:
		return

	#pairs=cmdline.split()
	for p in cmdline:
		tags=p.split(':')
		settings[tags[0]]=level_strings[tags[1]]


	for tr in tracers.values():
		if tr._full_name() in settings:
			tr.level=settings[tr._full_name()]
		elif not tr.parent:
			if 'all' in settings:
				tr.level = settings['all']
			else:
				tr.level = _tracer_info

class BaseTracer:
	def __init__(self, channel):
		self.channel = channel
		self.level = None

	def Warning( self, *msg ):
		self._start_output( msg, _tracer_info )

	def Info( self, *msg ):
		self._start_output( msg, _tracer_info )

	def Debug( self, *msg ):
		self._start_output( msg, _tracer_debug )

	def out( self, *msg ):
		self._start_output( msg, _tracer_debug )

	def Trace( self, *msg ):
		self._start_output( msg, _tracer_trace )

	def visible( self, level ):
		return self.level >= level or not self.level

	def _start_output( self, msg, level ):
		channels = []
		self._print_output( msg, level, channels )

	def _print_output( self, msg, level, channels, visible=None ):
		  if not visible: return
			#prepand your own channel
		  channels.insert(0,self.channel)
			prefix=".".join(channels)+": "
			astr=''
			for s in msg:
					try:
						astr+=" "+'%s'%str(s)
					except TypeError:
						#couldn't convert to string... maybe it is an iterable?
						try:
							print 'used tracer sub-code!'
							for sub in s:
								astr+=" "+'%s'%sub
						except:
							print 'TRACER: (cannot represent):',s
			for line in astr.split('\n'):
				print prefix, line

	def set_priority( self, level ):
		self.level = level

	def __repr__(self):
		return '(%s, %s)'%(self.channel,self.level)

class _Tracer( BaseTracer ):
	def __init__(self, subchannel, parent=None ):
		channels=subchannel.split('.')
		head='.'.join(channels[0:-1])
		tail=channels[-1]
		if head:
			parent=tracers.setdefault( head, _Tracer(head) )
		self.parent = parent
		BaseTracer.__init__( self, tail )

		#set levels from cmd-line settings
		if subchannel in settings:
			self.level=settings[ subchannel ]
		elif not self.parent:
			if 'all' in settings:
				self.level = settings['all']
			else:
				self.level = _tracer_info

	def _full_name(self):
		str=self.channel
		if self.parent:
			str=self.parent._full_name()+'.'+str
		return str

	def __repr__(self):
		return '(%s, %s, %s)'%(self.channel,self.parent, self.level)

	def _print_output( self, msg, level, channels, visible=None ):
		if visible==None and self.level:
			visible=self.visible(level)
		if self.parent:
			channels.insert(0,self.channel)
			self.parent._print_output( msg, level, channels, visible )
		else:
			BaseTracer._print_output( self, msg, level, channels, visible )

#User level Tracer class
class Tracer( _Tracer ):
	def __init__(self, subchannel, parent=None ):
		if parent:
			subchannel=parent._full_name()+'.'+subchannel
		self.my_tracer = tracers.setdefault( subchannel, _Tracer( subchannel ) )


	def _start_output( self, msg, level ):
		global warned_on_uninitialized
		channels = []
		if not initialized and not warned_on_uninitialized:
			print '. '*50
			print '[WARNING] Tracer system has not been initialized. Call library.init(args) at beginning of main program'
			print '. '*50
			warned_on_uninitialized=True

		self.my_tracer._print_output( msg, level, channels )

	def change_priority( self, level ):
		self.my_tracer.set_priority( level )

	def _full_name(self):
		return self.my_tracer._full_name()
