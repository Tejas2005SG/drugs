import React from 'react';
import {useAuthStore} from '../../Store/auth.store.js'

function VoicetoTextNotes() {
  const {user} = useAuthStore();
  console.log(user);

  return (
    <div>VoicetoTextNotes</div>
  )
}

export default VoicetoTextNotes