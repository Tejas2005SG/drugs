import jwt from "jsonwebtoken";

export const generateTokenAndCookie = (res, userId) => {
	const token = jwt.sign({ userId }, process.env.JWT_SECRET, {
		expiresIn: "7d",
	});

	res.cookie("token", token, {
		httpOnly: true,
		secure: true,
		sameSite: true,
		
	});

	return token;
};