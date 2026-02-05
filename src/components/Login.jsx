import React, { useState } from 'react';

const Login = ({ onLogin }) => {
    const [formData, setFormData] = useState({
        email: '',
        password: ''
    });

    const handleChange = (e) => {
        setFormData({
            ...formData,
            [e.target.name]: e.target.value
        });
    };

    const handleSubmit = (e) => {
        e.preventDefault();
        console.log('Login attempt:', formData);
        // In a real app, you would validate and make an API call here
        // For now, we simulate a successful login
        if (onLogin) onLogin(formData);
    };

    return (
        <div style={styles.container}>
            <div style={styles.card}>
                <div style={styles.header}>
                    <h2 style={styles.title}>Welcome Back</h2>
                    <p style={styles.subtitle}>Sign in to continue your research</p>
                </div>

                <form onSubmit={handleSubmit} style={styles.form}>
                    <div style={styles.inputGroup}>
                        <label style={styles.label}>Email Address</label>
                        <input
                            type="email"
                            name="email"
                            value={formData.email}
                            onChange={handleChange}
                            placeholder="Enter your email"
                            style={styles.input}
                            required
                        />
                    </div>

                    <div style={styles.inputGroup}>
                        <label style={styles.label}>Password</label>
                        <input
                            type="password"
                            name="password"
                            value={formData.password}
                            onChange={handleChange}
                            placeholder="Enter your password"
                            style={styles.input}
                            required
                        />
                    </div>

                    <button type="submit" style={styles.button}>
                        Sign In
                    </button>
                </form>
            </div>
        </div>
    );
};

const styles = {
    container: {
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        minHeight: '600px',
        padding: '20px',
        animation: 'fadeIn 0.5s ease-out',
    },
    card: {
        width: '100%',
        maxWidth: '400px',
        backgroundColor: '#ffffff',
        borderRadius: '16px',
        boxShadow: '0 4px 20px rgba(0, 0, 0, 0.05)',
        padding: '40px',
        border: '1px solid #f3f4f6',
    },
    header: {
        textAlign: 'center',
        marginBottom: '32px',
    },
    title: {
        fontSize: '28px',
        fontWeight: '700',
        color: '#111827',
        margin: '0 0 8px 0',
        fontFamily: "'DM Sans', sans-serif",
    },
    subtitle: {
        fontSize: '15px',
        color: '#6b7280',
        margin: 0,
        lineHeight: '1.5',
        fontFamily: "'DM Sans', sans-serif",
    },
    form: {
        display: 'flex',
        flexDirection: 'column',
        gap: '20px',
    },
    inputGroup: {
        display: 'flex',
        flexDirection: 'column',
        gap: '8px',
    },
    label: {
        fontSize: '14px',
        fontWeight: '600',
        color: '#374151',
        fontFamily: "'DM Sans', sans-serif",
    },
    input: {
        padding: '12px 16px',
        borderRadius: '8px',
        border: '1px solid #e5e7eb',
        fontSize: '15px',
        backgroundColor: '#f9fafb',
        transition: 'all 0.2s ease',
        outline: 'none',
        fontFamily: "'DM Sans', sans-serif",
    },
    button: {
        marginTop: '12px',
        padding: '14px',
        backgroundColor: '#16a34a',
        color: '#ffffff',
        border: 'none',
        borderRadius: '8px',
        fontSize: '16px',
        fontWeight: '600',
        cursor: 'pointer',
        transition: 'background-color 0.2s ease',
        fontFamily: "'DM Sans', sans-serif",
    }
};

export default Login;
