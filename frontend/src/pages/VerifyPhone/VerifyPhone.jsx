import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { useAuthStore } from '../../Store/auth.store.js';
import { toast } from 'react-hot-toast';

function VerifyPhone() {
  const [otp, setOtp] = useState('');
  const { verifyPhone, loading, phoneNumber } = useAuthStore();
  const navigate = useNavigate();

  const handleSubmit = async (e) => {
    e.preventDefault();

    // Validate OTP
    if (!otp || otp.length !== 6 || !/^\d+$/.test(otp)) {
      return toast.error('Please enter a valid 6-digit OTP');
    }

    if (!phoneNumber) {
      return toast.error('Phone number not found. Please sign up again.');
    }

    try {
      // Send both phoneNumber and otp to the backend
      await verifyPhone({ phoneNumber, otp });
      toast.success('Phone number verified successfully!');
      setOtp('');
      navigate('/dashboard'); // Navigate to dashboard on success
    } catch (error) {
      const errorMessage = error.response?.data?.message || 'OTP verification failed';
      toast.error(errorMessage);
    }
  };

  return (
    <div className="min-h-screen flex items-center justify-center bg-gray-100">
      <div className="bg-white p-8 rounded-lg shadow-md w-full max-w-md">
        <h2 className="text-2xl font-bold mb-6 text-center">Verify Phone Number</h2>
        <p className="text-center mb-4">
          An OTP has been sent to {phoneNumber || 'your phone number'}
        </p>
        <form onSubmit={handleSubmit}>
          <div className="mb-4">
            <label htmlFor="otp" className="block text-sm font-medium text-gray-700">
              Enter OTP
            </label>
            <input
              type="text"
              id="otp"
              value={otp}
              onChange={(e) => setOtp(e.target.value)}
              className="mt-1 block w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-indigo-500 focus:border-indigo-500"
              placeholder="Enter 6-digit OTP"
              maxLength={6}
              disabled={loading}
            />
          </div>
          <button
            type="submit"
            className="w-full bg-indigo-600 text-white py-2 px-4 rounded-md hover:bg-indigo-700 focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2 disabled:bg-indigo-400"
            disabled={loading}
          >
            {loading ? (
              <>
                <svg
                  className="animate-spin -ml-1 mr-2 h-5 w-5 text-white inline"
                  xmlns="http://www.w3.org/2000/svg"
                  fill="none"
                  viewBox="0 0 24 24"
                >
                  <circle
                    className="opacity-25"
                    cx="12"
                    cy="12"
                    r="10"
                    stroke="currentColor"
                    strokeWidth="4"
                  ></circle>
                  <path
                    className="opacity-75"
                    fill="currentColor"
                    d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                  ></path>
                </svg>
                Verifying...
              </>
            ) : (
              'Verify OTP'
            )}
          </button>
        </form>
      </div>
    </div>
  );
}

export default VerifyPhone;