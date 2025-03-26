"use client"

import { useState, useEffect } from "react";
import io from "socket.io-client";
import axios from "axios";
import { useAuthStore } from "../../Store/auth.store.js";

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:5000";

const axiosInstance = axios.create({
  baseURL: API_BASE_URL,
  withCredentials: true,
});

const socket = io(API_BASE_URL, {
  withCredentials: true,
  reconnection: true,
});

function Message() {
  const [messages, setMessages] = useState([]);
  const [newMessage, setNewMessage] = useState("");
  const [users, setUsers] = useState([]);
  const [selectedUser, setSelectedUser] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const { user, checkAuth, checkingAuth } = useAuthStore();

  useEffect(() => {
    const initializeChat = async () => {
      setLoading(true);
      await checkAuth();
      const currentUser = useAuthStore.getState().user;

      if (!currentUser) {
        setError("Authentication failed. Please log in.");
        setLoading(false);
        return;
      }

      const userId = currentUser._id;
      socket.emit("join", userId);
      console.log(`Client: User ${userId} joining room`);

      socket.on("connect", () => {
        console.log(`Client: Connected to socket server for user ${userId}`);
      });

      socket.on("newMessage", (message) => {
        console.log("Client: New message received:", message);
        setMessages((prev) => {
          const isRelevant =
            (message.sender._id === userId && message.receiver._id === selectedUser?._id) ||
            (message.sender._id === selectedUser?._id && message.receiver._id === userId);
          if (!isRelevant) {
            console.log("Client: Message ignored - not for current conversation");
            return prev;
          }
          if (prev.some((msg) => msg._id === message._id)) {
            console.log("Client: Duplicate message ignored");
            return prev;
          }
          console.log("Client: Adding message to state");
          return [...prev, message];
        });
      });

      socket.on("connect_error", (err) => {
        console.error("Client: Socket connection error:", err.message);
      });

      try {
        const response = await axiosInstance.get("/api/message/users/list");
        const validUsers = response.data.filter((u) => u.username && typeof u.username === "string");
        setUsers(validUsers);
      } catch (err) {
        setError(err.response?.data?.message || "Failed to fetch users");
      } finally {
        setLoading(false);
      }
    };

    initializeChat();

    return () => {
      socket.off("connect");
      socket.off("newMessage");
      socket.off("connect_error");
    };
  }, [checkAuth, selectedUser]);

  useEffect(() => {
    const fetchMessages = async () => {
      if (!selectedUser || !user) return;
      try {
        setLoading(true);
        const response = await axiosInstance.get(`/api/message/${selectedUser._id}`);
        setMessages(response.data);
      } catch (err) {
        setError(err.response?.data?.message || "Failed to fetch messages");
      } finally {
        setLoading(false);
      }
    };

    fetchMessages();
  }, [selectedUser, user]);

  const handleSendMessage = async (e) => {
    e.preventDefault();
    if (!newMessage.trim() || !selectedUser || !user) return;

    try {
      const response = await axiosInstance.post("/api/message/send", {
        receiverId: selectedUser._id,
        content: newMessage,
      });

      if (response.status === 201) {
        console.log("Client: Message sent successfully:", response.data);
        setNewMessage(""); // Clear input, let Socket.IO update messages
      }
    } catch (error) {
      setError(error.response?.data?.message || "Failed to send message");
    }
  };

  useEffect(() => {
    const messageContainer = document.getElementById("message-container");
    if (messageContainer) {
      messageContainer.scrollTop = messageContainer.scrollHeight;
    }
  }, [messages]);

  // JSX remains unchanged
  if (checkingAuth || (loading && !selectedUser)) {
    return (
      <div className="flex items-center justify-center h-screen bg-gray-50">
        <div className="text-center space-y-4">
          <div className="inline-block animate-spin rounded-full h-12 w-12 border-4 border-t-transparent border-blue-600"></div>
          <p className="text-gray-600 text-lg">
            {checkingAuth ? "Verifying authentication..." : "Loading messages..."}
          </p>
        </div>
      </div>
    );
  }

  if (!user) {
    return (
      <div className="flex items-center justify-center h-screen bg-gray-50">
        <div className="text-center p-8 bg-white rounded-2xl shadow-lg max-w-md space-y-6">
          <div className="space-y-4">
            <div className="text-blue-600 mx-auto animate-pulse">
              <svg className="w-20 h-20" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth="1.5"
                  d="M12 15v2m-6 4h12a2 2 0 002-2v-6a2 2 0 00-2-2H6a2 2 0 00-2 2v6a2 2 0 002 2zm10-10V7a4 4 0 00-8 0v4h8z"
                />
              </svg>
            </div>
            <h2 className="text-3xl font-bold text-gray-800">Access Restricted</h2>
            <p className="text-gray-600 text-lg">Sign in to use messaging</p>
          </div>
          <button
            className="w-full px-6 py-3 bg-blue-600 text-white rounded-xl 
                      hover:bg-blue-700 transition-all duration-300 font-semibold 
                      shadow-md hover:shadow-lg"
            onClick={() => (window.location.href = "/login")}
          >
            Continue to Login
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="flex h-screen max-h-screen overflow-hidden bg-gray-100">
      {/* Sidebar */}
      <div className="w-1/4 min-w-[250px] bg-white border-r border-gray-200 flex flex-col">
        <div className="p-4 border-b border-gray-200 bg-white sticky top-0 z-10">
          <h2 className="text-xl font-semibold text-gray-800">Messages</h2>
        </div>
        <div className="flex-1 overflow-y-auto scrollbar-thin scrollbar-thumb-gray-300 scrollbar-track-transparent">
          {users.length === 0 ? (
            <div className="flex flex-col items-center justify-center h-full p-6 text-center">
              <div className="w-16 h-16 bg-gray-100 rounded-full flex items-center justify-center mb-4">
                <svg className="w-8 h-8 text-gray-400" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path
                    strokeLinecap="round"
                    strokeLinejoin="round"
                    strokeWidth="2"
                    d="M16 7a4 4 0 11-8 0 4 4 0 018 0zM12 14a7 7 0 00-7 7h14a7 7 0 00-7-7z"
                  />
                </svg>
              </div>
              <p className="text-gray-500">No users available</p>
              <p className="text-sm text-gray-400 mt-2">Users will appear here once they join</p>
            </div>
          ) : (
            users.map((u) => (
              <div
                key={u._id}
                onClick={() => setSelectedUser(u)}
                className={`flex items-center p-4 cursor-pointer hover:bg-gray-50 transition-colors ${
                  selectedUser?._id === u._id
                    ? "bg-blue-50 border-l-4 border-blue-500"
                    : "border-l-4 border-transparent"
                }`}
              >
                <div className="w-12 h-12 rounded-full bg-gradient-to-r from-blue-400 to-indigo-500 flex items-center justify-center text-white font-medium shadow-sm">
                  {u.profilePicture ? (
                    <img
                      src={u.profilePicture || "/placeholder.svg"}
                      alt={u.username}
                      className="w-full h-full rounded-full object-cover"
                    />
                  ) : (
                    <span>{u.username ? u.username[0]?.toUpperCase() : "?"}</span>
                  )}
                </div>
                <div className="ml-3 flex-1 overflow-hidden">
                  <p className="text-gray-800 font-medium truncate">{u.username || "Unknown User"}</p>
                  <p className="text-gray-500 text-sm truncate">Click to start chatting</p>
                </div>
                <div className="w-2 h-2 rounded-full bg-green-500 ml-2"></div>
              </div>
            ))
          )}
        </div>
      </div>

      {/* Chat Window */}
      <div className="flex-1 flex flex-col h-full">
        {selectedUser ? (
          <>
            {/* Chat Header */}
            <div className="p-4 bg-white border-b border-gray-200 flex items-center sticky top-0 z-10 shadow-sm">
              <div className="w-10 h-10 rounded-full bg-gradient-to-r from-blue-400 to-indigo-500 flex items-center justify-center text-white font-medium">
                {selectedUser.profilePicture ? (
                  <img
                    src={selectedUser.profilePicture || "/placeholder.svg"}
                    alt={selectedUser.username}
                    className="w-full h-full rounded-full object-cover"
                  />
                ) : (
                  <span>{selectedUser.username ? selectedUser.username[0]?.toUpperCase() : "?"}</span>
                )}
              </div>
              <div className="ml-3">
                <h2 className="text-lg font-semibold text-gray-800">{selectedUser.username || "Unknown User"}</h2>
                <p className="text-xs text-gray-500">Online</p>
              </div>
            </div>

            {/* Messages */}
            <div
              id="message-container"
              className="flex-1 p-4 overflow-y-auto bg-gray-50 scrollbar-thin scrollbar-thumb-gray-300 scrollbar-track-transparent"
              style={{ maxHeight: "calc(100vh - 140px)" }}
            >
              {loading && !messages.length ? (
                <div className="space-y-4">
                  {[...Array(3)].map((_, i) => (
                    <div key={i} className="animate-pulse flex space-x-4">
                      <div className="h-10 w-10 bg-gray-200 rounded-full"></div>
                      <div className="flex-1 space-y-2">
                        <div className="h-4 bg-gray-200 rounded w-1/4"></div>
                        <div className="h-4 bg-gray-200 rounded w-3/4"></div>
                      </div>
                    </div>
                  ))}
                </div>
              ) : error ? (
                <div className="p-4 bg-red-50 rounded-xl flex items-center gap-4 max-w-md mx-auto my-4">
                  <svg
                    className="w-8 h-8 text-red-600 flex-shrink-0"
                    fill="none"
                    stroke="currentColor"
                    viewBox="0 0 24 24"
                  >
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth="2"
                      d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-3L13.732 4c-.77-1.333-2.694-1.333-3.464 0L3.34 16c-.77 1.333.192 3 1.732 3z"
                    />
                  </svg>
                  <div>
                    <h3 className="font-medium text-red-700">Error</h3>
                    <p className="text-sm text-red-600 mt-1">{error}</p>
                  </div>
                </div>
              ) : messages.length === 0 ? (
                <div className="flex flex-col items-center justify-center h-full">
                  <div className="w-24 h-24 bg-blue-50 rounded-full flex items-center justify-center mb-4">
                    <svg className="w-12 h-12 text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path
                        strokeLinecap="round"
                        strokeLinejoin="round"
                        strokeWidth="1.5"
                        d="M8 12h.01M12 12h.01M16 12h.01M21 12c0 4.418-4.03 8-9 8a9.863 9.863 0 01-4.255-.949L3 20l1.395-3.72C3.512 15.042 3 13.574 3 12c0-4.418 4.03-8 9-8s9 3.582 9 8z"
                      />
                    </svg>
                  </div>
                  <p className="text-gray-600 font-medium">
                    Start a conversation with {selectedUser.username || "this user"}!
                  </p>
                  <p className="text-gray-500 text-sm mt-2">Your messages will appear here</p>
                </div>
              ) : (
                <div className="space-y-4">
                  {messages.map((msg, index) => {
                    const isCurrentUser = msg.sender._id === user._id;
                    const showAvatar = index === 0 || messages[index - 1].sender._id !== msg.sender._id;

                    return (
                      <div key={msg._id} className={`flex ${isCurrentUser ? "justify-end" : "justify-start"}`}>
                        {!isCurrentUser && showAvatar && (
                          <div className="w-8 h-8 rounded-full bg-gradient-to-r from-purple-400 to-pink-500 flex items-center justify-center text-white text-xs mr-2 flex-shrink-0">
                            {selectedUser.username ? selectedUser.username[0]?.toUpperCase() : "?"}
                          </div>
                        )}

                        <div
                          className={`max-w-xs md:max-w-md lg:max-w-lg p-3 rounded-2xl ${
                            isCurrentUser
                              ? "bg-gradient-to-r from-blue-500 to-blue-600 text-white rounded-tr-none shadow-md"
                              : "bg-white text-gray-800 rounded-tl-none shadow-sm border border-gray-100"
                          }`}
                        >
                          <p className="whitespace-pre-wrap break-words">{msg.content}</p>
                          <p className={`text-xs mt-1 text-right ${isCurrentUser ? "text-blue-100" : "text-gray-500"}`}>
                            {new Date(msg.timestamp).toLocaleTimeString([], {
                              hour: "numeric",
                              minute: "2-digit",
                              hour12: true,
                            })}
                          </p>
                        </div>

                        {isCurrentUser && showAvatar && (
                          <div className="w-8 h-8 rounded-full bg-gradient-to-r from-blue-400 to-indigo-500 flex items-center justify-center text-white text-xs ml-2 flex-shrink-0">
                            {user.username ? user.username[0]?.toUpperCase() : "?"}
                          </div>
                        )}
                      </div>
                    );
                  })}
                </div>
              )}
            </div>

            {/* Message Input */}
            <div className="p-4 bg-white border-t border-gray-200 sticky bottom-0">
              <form onSubmit={handleSendMessage} className="flex items-center gap-2">
                <input
                  type="text"
                  value={newMessage}
                  onChange={(e) => setNewMessage(e.target.value)}
                  placeholder="Type a message..."
                  className="flex-1 px-4 py-3 border border-gray-300 rounded-full focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
                />
                <button
                  type="submit"
                  className="p-3 bg-gradient-to-r from-blue-500 to-blue-600 text-white rounded-full hover:from-blue-600 hover:to-blue-700 transition-all disabled:opacity-50 disabled:cursor-not-allowed shadow-md hover:shadow-lg"
                  disabled={!newMessage.trim()}
                >
                  <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth="2"
                      d="M12 19l9 2-9-18-9 18 9-2zm0 0v-8"
                    />
                  </svg>
                </button>
              </form>
            </div>
          </>
        ) : (
          <div className="flex flex-col items-center justify-center h-full bg-gray-50 p-6">
            <div className="w-24 h-24 bg-blue-50 rounded-full flex items-center justify-center mb-6">
              <svg className="w-12 h-12 text-blue-500" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth="1.5"
                  d="M8 10h.01M12 10h.01M16 10h.01M9 16H5a2 2 0 01-2-2V6a2 2 0 012-2h14a2 2 0 012 2v8a2 2 0 01-2 2h-5l-5 5v-5z"
                />
              </svg>
            </div>
            <h3 className="text-xl font-semibold text-gray-800 mb-2">Your Messages</h3>
            <p className="text-gray-500 text-center max-w-md">Select a user from the sidebar to start chatting</p>
          </div>
        )}
      </div>
    </div>
  );
}

export default Message;