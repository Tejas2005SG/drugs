import React from "react";

export const Spinner = () => {
  return (
    <div className="flex justify-center">
      <div className="animate-spin h-8 w-8 rounded-full border-4 border-primary-200 border-t-primary-500"></div>
    </div>
  );
};